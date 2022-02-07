# Module containing functions to read from the SQLite

from condenser import NumpyQuery
from contextlib import contextmanager
from enum import Enum
from functools import lru_cache
from pyproj import Transformer
from pyproj.crs import CRS
from sqlalchemy import cast
from sqlalchemy import inspect
from sqlalchemy import Integer
from sqlalchemy.orm import Session
from threedi_modelchecker.schema import ModelSchema
from threedi_modelchecker.threedi_database import ThreediDatabase
from threedi_modelchecker.threedi_model import models
from threedi_modelchecker.threedi_model.custom_types import IntegerEnum
from threedigrid_builder.base import GridSettings
from threedigrid_builder.base import Levees
from threedigrid_builder.base import Pumps
from threedigrid_builder.base import TablesSettings
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import InitializationType
from threedigrid_builder.constants import LineType
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import BoundaryConditions1D
from threedigrid_builder.grid import BoundaryConditions2D
from threedigrid_builder.grid import Channels
from threedigrid_builder.grid import ConnectedPoints
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import CrossSectionDefinitions
from threedigrid_builder.grid import CrossSectionLocations
from threedigrid_builder.grid import Culverts
from threedigrid_builder.grid import DemAverageAreas
from threedigrid_builder.grid import GridRefinements
from threedigrid_builder.grid import ImperviousSurfaces
from threedigrid_builder.grid import Obstacles
from threedigrid_builder.grid import Orifices
from threedigrid_builder.grid import Pipes
from threedigrid_builder.grid import Surfaces
from threedigrid_builder.grid import Weirs
from typing import Callable
from typing import ContextManager
from typing import Tuple

import numpy as np
import pathlib
import pygeos


__all__ = ["SQLite"]

# hardcoded source projection
SOURCE_EPSG = 4326

# before version 173, we get an error because "interception_file" is missing
MIN_SQLITE_VERSION = 173

# put some global defaults on datatypes
NumpyQuery.default_numpy_settings[Integer] = {"dtype": np.int32, "null": -9999}
NumpyQuery.default_numpy_settings[IntegerEnum] = {
    **NumpyQuery.default_numpy_settings[Integer],
    "sql_cast": lambda x: cast(x, Integer),
}

# To parse CalculationPoint.user_ref
CALCULATION_POINT_CONTENT_TYPE_MAP = {
    "v2_channel": ContentType.TYPE_V2_CHANNEL,
    "v2_pipe": ContentType.TYPE_V2_PIPE,
    "v2_culvert": ContentType.TYPE_V2_CULVERT,
    "v2_manhole": ContentType.TYPE_V2_MANHOLE,
    "v2_1d_boundary_conditions": ContentType.TYPE_V2_1D_BOUNDARY_CONDITIONS,
}


def parse_connected_point_user_ref(user_ref: str) -> Tuple[ContentType, int, int]:
    """Return content_type, content_id, node_number from a user_ref.

    Raises Exception for various parse errors.

    Example
    -------
    >>> parse_connected_point_user_ref("201#123#v2_channels#4)
    ContentType.TYPE_V2_CHANNEL, 123, 4
    """
    _, id_str, type_str, node_str = user_ref.split("#")
    return CALCULATION_POINT_CONTENT_TYPE_MAP[type_str], int(id_str), int(node_str)


def _object_as_dict(obj) -> dict:
    """Convert SQLAlchemy object to dict, casting Enums"""
    result = {}
    if obj is None:
        return result
    for c in inspect(obj).mapper.column_attrs:
        val = getattr(obj, c.key)
        if isinstance(val, Enum):
            val = val.value
        result[c.key] = val
    return result


def _set_initialization_type(
    dct, global_field, file_field=None, type_field=None, default=None
):
    """Set the InitializationType depending on global_field and file_field."""
    if not file_field:
        file_field = f"{global_field}_file"
    if not type_field:
        type_field = f"{global_field}_type"

    # If the ``file_field`` contains a value, the initialization type will be changed to
    # the ``default_type``, if supplied.
    if dct[file_field]:
        if default is not None:
            dct[type_field] = default
    # If there is no file, check if the global value is not None
    elif dct[global_field] is not None:
        dct[type_field] = InitializationType.GLOBAL
    else:
        # No file, no global value
        dct[type_field] = None


class SQLite:
    def __init__(self, path: pathlib.Path):
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")
        path = str(path)
        sqlite_settings = {"db_path": path, "db_file": path}
        self.db = ThreediDatabase(
            connection_settings=sqlite_settings, db_type="spatialite"
        )
        self._epsg_code = None  # for reproject()

        version = self.get_version()
        if version < MIN_SQLITE_VERSION:
            raise SchematisationError(f"Too old sqlite version {version}.")

    def get_version(self) -> int:
        # check version
        schema = ModelSchema(self.db)
        return schema.get_version()

    @contextmanager
    def get_session(self) -> ContextManager[Session]:
        """A context manager that yields an SQLAlchemy session.

        The session is closed af the context manager exit. No commit or rollback
        is done. It is meant for read-only access, but writing is not prohibited.

        Returns:
          SQLAlchemy.orm.Session
        """
        session = self.db.get_session(query_cls=NumpyQuery)
        try:
            yield session
        finally:
            session.close()

    def get_settings(self) -> dict:
        """Return the settings relevant for makegrid and maketables.

        Returns:
           dict with epsg_code, model_name, grid_settings, make_table_settings
        """

        with self.get_session() as session:
            global_ = session.query(models.GlobalSetting).order_by("id").first()
            if global_.groundwater_settings_id is not None:
                groundwater = _object_as_dict(
                    session.query(models.GroundWater)
                    .filter_by(id=global_.groundwater_settings_id)
                    .one()
                )
            else:
                groundwater = {}
            if global_.interflow_settings_id is not None:
                interflow = _object_as_dict(
                    session.query(models.Interflow)
                    .filter_by(id=global_.interflow_settings_id)
                    .one()
                )
            else:
                interflow = {}
            if global_.simple_infiltration_settings_id is not None:
                infiltration = _object_as_dict(
                    session.query(models.SimpleInfiltration)
                    .filter_by(id=global_.simple_infiltration_settings_id)
                    .one()
                )
            else:
                infiltration = {}
            global_ = _object_as_dict(global_)

        # record if there is a DEM file to be expected
        # Note: use_2d_flow only determines whether there are flow lines
        global_["use_2d"] = bool(global_["dem_file"])

        # set/adapt initialization types to include information about file presence
        NO_AGG = InitializationType.NO_AGG
        AVERAGE = InitializationType.AVERAGE
        _set_initialization_type(
            global_, "frict_coef", default=AVERAGE if global_["frict_avg"] else NO_AGG
        )
        _set_initialization_type(
            global_,
            "interception_global",
            file_field="interception_file",
            type_field="interception_type",
            default=NO_AGG,
        )
        if interflow:
            _set_initialization_type(interflow, "porosity", default=NO_AGG)
            _set_initialization_type(
                interflow, "hydraulic_conductivity", default=AVERAGE
            )
        if infiltration:
            _set_initialization_type(infiltration, "infiltration_rate", default=NO_AGG)
            # max_infiltration_capacity_file has no corresponding global value!
            infiltration["max_infiltration_capacity_file"] = infiltration.get(
                "max_infiltration_capacity_file"
            )
            if (
                infiltration["max_infiltration_capacity_file"] is not None
                and infiltration["max_infiltration_capacity_file"] != ""
            ):
                infiltration["max_infiltration_capacity_type"] = NO_AGG
            else:
                infiltration["max_infiltration_capacity_type"] = None

        if groundwater:
            # default is what the user supplied (MIN/MAX/AVERAGE)
            _set_initialization_type(groundwater, "groundwater_impervious_layer_level")
            _set_initialization_type(groundwater, "phreatic_storage_capacity")
            _set_initialization_type(groundwater, "equilibrium_infiltration_rate")
            _set_initialization_type(groundwater, "initial_infiltration_rate")
            _set_initialization_type(groundwater, "infiltration_decay_period")
            _set_initialization_type(groundwater, "groundwater_hydro_connectivity")

        grid_settings = GridSettings.from_dict(global_)
        tables_settings = TablesSettings.from_dict(
            {**groundwater, **interflow, **infiltration, **global_}
        )
        return {
            "epsg_code": global_["epsg_code"],
            "model_name": global_["name"],
            "grid_settings": grid_settings,
            "tables_settings": tables_settings,
        }

    @property
    def epsg_code(self) -> int:
        if self._epsg_code is None:
            self._epsg_code = self.get_settings()["epsg_code"]
        return self._epsg_code

    @epsg_code.setter
    def epsg_code(self, value: int):
        self._epsg_code = value

    def reproject(self, geometries: np.ndarray) -> np.ndarray:
        """Reproject geometries from 4326 to the EPSG in the settings.

        Notes:
          pygeos+pyproj is approx 2x faster than spatialite

        Args:
          geometries (ndarray of pygeos.Geometry): geometries in EPSG 4326
        """
        target_epsg = self.epsg_code
        func = _get_reproject_func(SOURCE_EPSG, target_epsg)
        return pygeos.apply(geometries, func)

    def get_surfaces(self) -> Surfaces:
        with self.get_session() as session:
            arr = (
                session.query(
                    models.Surface.id.label("surface_id"),
                    models.Surface.function,
                    models.Surface.code,
                    models.Surface.display_name,
                    models.Surface.nr_of_inhabitants,
                    models.Surface.area,
                    models.Surface.dry_weather_flow,
                    models.Surface.the_geom,
                    models.SurfaceParameter.outflow_delay,
                    models.SurfaceParameter.surface_layer_thickness,
                    models.SurfaceParameter.infiltration,
                    models.SurfaceParameter.max_infiltration_capacity,
                    models.SurfaceParameter.min_infiltration_capacity,
                    models.SurfaceParameter.infiltration_decay_constant,
                    models.SurfaceParameter.infiltration_recovery_constant,
                    models.ConnectionNode.id.label("connection_node_id"),
                    models.ConnectionNode.the_geom.label("connection_node_the_geom"),
                    models.SurfaceMap.percentage,
                )
                .select_from(models.Surface)
                .join(models.SurfaceParameter)
                .join(
                    models.SurfaceMap, models.SurfaceMap.surface_id == models.Surface.id
                )
                .join(
                    models.ConnectionNode,
                    models.SurfaceMap.connection_node_id == models.ConnectionNode.id,
                )
                .order_by(models.Surface.id)
                .as_structarray()
            )

        # reproject
        arr["the_geom"] = self.reproject(arr["the_geom"])
        arr["connection_node_the_geom"] = self.reproject(
            arr["connection_node_the_geom"]
        )

        return Surfaces(
            id=np.arange(0, len(arr["surface_id"] + 1), dtype=int),
            **{name: arr[name] for name in arr.dtype.names},
        )

    def get_impervious_surfaces(self) -> ImperviousSurfaces:
        with self.get_session() as session:
            arr = (
                session.query(
                    models.ImperviousSurface.id.label("surface_id"),
                    models.ImperviousSurface.code,
                    models.ImperviousSurface.display_name,
                    models.ImperviousSurface.surface_inclination,
                    models.ImperviousSurface.surface_class,
                    models.ImperviousSurface.surface_sub_class,
                    models.ImperviousSurface.nr_of_inhabitants,
                    models.ImperviousSurface.area,
                    models.ImperviousSurface.dry_weather_flow,
                    models.ImperviousSurface.the_geom,
                    models.ConnectionNode.id.label("connection_node_id"),
                    models.ConnectionNode.the_geom.label("connection_node_the_geom"),
                    models.ImperviousSurfaceMap.percentage,
                )
                .select_from(models.ImperviousSurface)
                .join(models.ImperviousSurfaceMap)
                .join(
                    models.ConnectionNode,
                    models.ImperviousSurfaceMap.connection_node_id
                    == models.ConnectionNode.id,
                )
                .order_by(models.ImperviousSurface.id)
                .as_structarray()
            )

        # convert enums to values
        arr["surface_class"] = [x.value for x in arr["surface_class"]]
        arr["surface_inclination"] = [x.value for x in arr["surface_inclination"]]

        # reproject
        arr["the_geom"] = self.reproject(arr["the_geom"])
        arr["connection_node_the_geom"] = self.reproject(
            arr["connection_node_the_geom"]
        )

        return ImperviousSurfaces(
            id=np.arange(0, len(arr["surface_id"] + 1), dtype=int),
            **{name: arr[name] for name in arr.dtype.names},
        )

    def get_boundary_conditions_1d(self) -> BoundaryConditions1D:
        """Return BoundaryConditions1D"""
        with self.get_session() as session:
            arr = (
                session.query(
                    models.BoundaryCondition1D.id,
                    models.BoundaryCondition1D.boundary_type,
                    models.BoundaryCondition1D.connection_node_id,
                )
                .order_by(models.BoundaryCondition1D.id)
                .as_structarray()
            )

        # transform to a BoundaryConditions1D object
        return BoundaryConditions1D(**{name: arr[name] for name in arr.dtype.names})

    def get_boundary_conditions_2d(self) -> BoundaryConditions2D:
        """Return BoundaryConditions2D"""
        with self.get_session() as session:
            arr = (
                session.query(
                    models.BoundaryConditions2D.id,
                    models.BoundaryConditions2D.boundary_type,
                    models.BoundaryConditions2D.the_geom,
                )
                .order_by(models.BoundaryConditions2D.id)
                .as_structarray()
            )
        arr["the_geom"] = self.reproject(arr["the_geom"])

        # transform to a BoundaryConditions1D object
        return BoundaryConditions2D(**{name: arr[name] for name in arr.dtype.names})

    def get_connected_points(self) -> ConnectedPoints:
        """Load connected points, join them directly with calculation points.

        Automatically ignores calculation points without connected points.
        """
        with self.get_session() as session:
            arr = (
                session.query(
                    models.ConnectedPoint.the_geom,
                    models.ConnectedPoint.id,
                    models.ConnectedPoint.exchange_level,
                    models.ConnectedPoint.levee_id,
                    models.CalculationPoint.user_ref,
                    models.CalculationPoint.id.label("calculation_point_id"),
                    models.Levee.crest_level,
                )
                .join(models.CalculationPoint)
                .outerjoin(models.Levee)
                .order_by(models.ConnectedPoint.id)
                .as_structarray()
            )

        # reproject
        arr["the_geom"] = self.reproject(arr["the_geom"])

        # replace -9999.0 with NaN in exchange_level
        arr["exchange_level"][arr["exchange_level"] == -9999.0] = np.nan

        # convert to columnar dict
        dct = {name: arr[name] for name in arr.dtype.names}

        # overwrite exchange_level with crest_level where no exchange_level supplied
        mask = np.isnan(arr["exchange_level"])
        dct["exchange_level"][mask] = dct.pop("crest_level")[mask]

        # parse user_ref
        content_type = np.empty_like(dct["id"])
        content_pk = np.empty_like(dct["id"])
        node_number = np.empty_like(dct["id"])
        for i, user_ref in enumerate(dct.pop("user_ref")):
            try:
                parsed = parse_connected_point_user_ref(user_ref)
            except Exception:
                raise SchematisationError(
                    f'Invalid user_ref in connected point {dct["id"][i]}: "{user_ref}".'
                )
            content_type[i], content_pk[i], node_number[i] = parsed

        # transform to a ConnectedPoints object
        return ConnectedPoints(
            content_type=content_type,
            content_pk=content_pk,
            node_number=node_number,
            **dct,
        )

    def get_levees(self) -> Levees:
        with self.get_session() as session:
            arr = (
                session.query(
                    models.Levee.id,
                    models.Levee.the_geom,
                    models.Levee.material,
                    models.Levee.crest_level,
                    models.Levee.max_breach_depth,
                )
                .order_by(models.Levee.id)
                .as_structarray()
            )

        # reproject
        arr["the_geom"] = self.reproject(arr["the_geom"])

        # replace -9999.0 with NaN in crest_level and max_breach_depth
        arr["crest_level"][arr["crest_level"] == -9999.0] = np.nan
        arr["max_breach_depth"][arr["max_breach_depth"] == -9999.0] = np.nan

        # transform to a Levees object
        return Levees(**{name: arr[name] for name in arr.dtype.names})

    def get_channels(self) -> Channels:
        """Return Channels"""
        with self.get_session() as session:
            arr = (
                session.query(
                    models.Channel.the_geom,
                    models.Channel.dist_calc_points,
                    models.Channel.id,
                    models.Channel.code,
                    models.Channel.connection_node_start_id,
                    models.Channel.connection_node_end_id,
                    models.Channel.calculation_type,
                    models.Channel.display_name,
                    models.Channel.zoom_category,
                )
                .order_by(models.Channel.id)
                .as_structarray()
            )

        arr["the_geom"] = self.reproject(arr["the_geom"])
        # map "old" calculation types (100, 101, 102, 105) to (0, 1, 2, 5)
        arr["calculation_type"][arr["calculation_type"] >= 100] -= 100

        # transform to a Channels object
        return Channels(**{name: arr[name] for name in arr.dtype.names})

    def get_connection_nodes(self) -> ConnectionNodes:
        """Return ConnectionNodes (which are enriched using the manhole table)"""
        with self.get_session() as session:
            arr = (
                session.query(
                    models.ConnectionNode.the_geom,
                    models.ConnectionNode.id,
                    models.ConnectionNode.code,
                    models.ConnectionNode.storage_area,
                    models.ConnectionNode.initial_waterlevel,
                    models.Manhole.id.label("manhole_id"),
                    models.Manhole.calculation_type,
                    models.Manhole.bottom_level,
                    models.Manhole.drain_level,
                    models.Manhole.manhole_indicator,
                    models.Manhole.surface_level,
                    models.Manhole.shape,
                    models.Manhole.width,
                    models.Manhole.display_name,
                    models.Manhole.zoom_category,
                )
                .join(models.ConnectionNode.manholes, isouter=True)
                .distinct(models.ConnectionNode.id)
                .order_by(models.ConnectionNode.id)
                .as_structarray()
            )

        arr["the_geom"] = self.reproject(arr["the_geom"])

        # replace -9999.0 with NaN in initial_waterlevel
        arr["initial_waterlevel"][arr["initial_waterlevel"] == -9999.0] = np.nan

        return ConnectionNodes(**{name: arr[name] for name in arr.dtype.names})

    def get_cross_section_definitions(self) -> CrossSectionDefinitions:
        """Return CrossSectionDefinitions"""
        with self.get_session() as session:
            arr = (
                session.query(
                    models.CrossSectionDefinition.id,
                    models.CrossSectionDefinition.code,
                    models.CrossSectionDefinition.shape,
                    models.CrossSectionDefinition.width,
                    models.CrossSectionDefinition.height,
                )
                .order_by(models.CrossSectionDefinition.id)
                .as_structarray()
            )

        # map shape 10 to 1 (circle) to match CrossSectionShape enum
        arr["shape"][arr["shape"] == 10] = 1
        # map shape 11 to 5 (tabulated rectangle) to match CrossSectionShape enum
        arr["shape"][arr["shape"] == 11] = 5
        # map shape 12 to 6 (tabulated trapezium) to match CrossSectionShape enum
        arr["shape"][arr["shape"] == 12] = 6

        # transform to a CrossSectionDefinitions object
        return CrossSectionDefinitions(**{name: arr[name] for name in arr.dtype.names})

    def get_cross_section_locations(self) -> CrossSectionLocations:
        """Return CrossSectionLocations"""
        with self.get_session() as session:
            arr = (
                session.query(
                    models.CrossSectionLocation.id,
                    models.CrossSectionLocation.code,
                    models.CrossSectionLocation.the_geom,
                    models.CrossSectionLocation.definition_id,
                    models.CrossSectionLocation.channel_id,
                    models.CrossSectionLocation.reference_level,
                    models.CrossSectionLocation.bank_level,
                    models.CrossSectionLocation.friction_type,
                    models.CrossSectionLocation.friction_value,
                )
                .order_by(models.CrossSectionLocation.id)
                .as_structarray()
            )

        arr["the_geom"] = self.reproject(arr["the_geom"])

        # transform to a CrossSectionLocations object
        return CrossSectionLocations(**{name: arr[name] for name in arr.dtype.names})

    def get_culverts(self) -> Culverts:
        """Return Culverts"""
        with self.get_session() as session:
            arr = (
                session.query(
                    models.Culvert.id,
                    models.Culvert.code,
                    models.Culvert.the_geom,
                    models.Culvert.dist_calc_points,
                    models.Culvert.connection_node_start_id,
                    models.Culvert.connection_node_end_id,
                    models.Culvert.calculation_type,
                    models.Culvert.cross_section_definition_id,
                    models.Culvert.invert_level_start_point,
                    models.Culvert.invert_level_end_point,
                    models.Culvert.discharge_coefficient_negative,
                    models.Culvert.discharge_coefficient_positive,
                    models.Culvert.friction_type,
                    models.Culvert.friction_value,
                    models.Culvert.display_name,
                    models.Culvert.zoom_category,
                )
                .order_by(models.Culvert.id)
                .as_structarray()
            )

        arr["the_geom"] = self.reproject(arr["the_geom"])

        # map friction_type 4 to friction_type 2 to match crosssectionlocation enum
        arr["friction_type"][arr["friction_type"] == 4] = 2

        # When no calculation type is provides we default to isolated
        arr["calculation_type"][
            arr["calculation_type"] == -9999
        ] = LineType.LINE_1D_ISOLATED
        # map "old" calculation types (100, 101, 102, 105) to (0, 1, 2, 5)
        arr["calculation_type"][arr["calculation_type"] >= 100] -= 100

        # transform to a CrossSectionLocations object
        return Culverts(**{name: arr[name] for name in arr.dtype.names})

    def get_grid_refinements(self) -> GridRefinements:
        """Return Gridrefinement and GridRefinementArea concatenated into one array."""
        with self.get_session() as session:
            arr1 = (
                session.query(
                    models.GridRefinement.the_geom,
                    models.GridRefinement.id,
                    models.GridRefinement.code,
                    models.GridRefinement.display_name,
                    models.GridRefinement.refinement_level,
                )
                .filter(
                    models.GridRefinement.the_geom.isnot(None),
                    models.GridRefinement.refinement_level.isnot(None),
                )
                .order_by(models.GridRefinement.id)
                .as_structarray()
            )
            arr2 = (
                session.query(
                    models.GridRefinementArea.the_geom,
                    models.GridRefinementArea.id,
                    models.GridRefinementArea.code,
                    models.GridRefinementArea.display_name,
                    models.GridRefinementArea.refinement_level,
                )
                .filter(
                    models.GridRefinementArea.the_geom.isnot(None),
                    models.GridRefinementArea.refinement_level.isnot(None),
                )
                .order_by(models.GridRefinementArea.id)
                .as_structarray()
            )
            arr = np.concatenate((arr1, arr2))

        # reproject
        arr["the_geom"] = self.reproject(arr["the_geom"])
        arr["id"] = np.arange(len(arr["refinement_level"]))

        return GridRefinements(**{name: arr[name] for name in arr.dtype.names})

    def get_dem_average_areas(self) -> DemAverageAreas:
        """Return DemAverageAreas"""
        with self.get_session() as session:
            arr = (
                session.query(
                    models.DemAverageArea.id,
                    models.DemAverageArea.the_geom,
                )
                .order_by(models.DemAverageArea.id)
                .as_structarray()
            )
            arr["the_geom"] = self.reproject(arr["the_geom"])

        return DemAverageAreas(**{name: arr[name] for name in arr.dtype.names})

    def get_obstacles(self) -> Obstacles:
        with self.get_session() as session:
            arr = (
                session.query(
                    models.Obstacle.the_geom,
                    models.Obstacle.id,
                    models.Obstacle.crest_level,
                )
                .order_by(models.Obstacle.id)
                .as_structarray()
            )

        # reproject
        arr["the_geom"] = self.reproject(arr["the_geom"])

        return Obstacles(**{name: arr[name] for name in arr.dtype.names})

    def get_orifices(self) -> Orifices:
        """Return Orifices"""
        with self.get_session() as session:
            arr = (
                session.query(
                    models.Orifice.id,
                    models.Orifice.code,
                    models.Orifice.connection_node_start_id,
                    models.Orifice.connection_node_end_id,
                    models.Orifice.crest_level,
                    models.Orifice.crest_type,
                    models.Orifice.cross_section_definition_id,
                    models.Orifice.discharge_coefficient_negative,
                    models.Orifice.discharge_coefficient_positive,
                    models.Orifice.friction_type,
                    models.Orifice.friction_value,
                    models.Orifice.display_name,
                    models.Orifice.zoom_category,
                    models.Orifice.sewerage,
                )
                .order_by(models.Orifice.id)
                .as_structarray()
            )

        # map friction_type 4 to friction_type 2 to match crosssectionlocation enum
        arr["friction_type"][arr["friction_type"] == 4] = 2

        return Orifices(**{name: arr[name] for name in arr.dtype.names})

    def get_pipes(self) -> Pipes:
        """Return Pipes"""
        with self.get_session() as session:
            arr = (
                session.query(
                    models.Pipe.id,
                    models.Pipe.code,
                    models.Pipe.sewerage_type,
                    models.Pipe.calculation_type,
                    models.Pipe.invert_level_start_point,
                    models.Pipe.invert_level_end_point,
                    models.Pipe.friction_type,
                    models.Pipe.friction_value,
                    models.Pipe.dist_calc_points,
                    models.Pipe.connection_node_start_id,
                    models.Pipe.connection_node_end_id,
                    models.Pipe.cross_section_definition_id,
                    models.Pipe.display_name,
                    models.Pipe.zoom_category,
                    models.Pipe.material,
                )
                .order_by(models.Pipe.id)
                .as_structarray()
            )

        # map friction_type 4 to friction_type 2 to match crosssectionlocation enum
        arr["friction_type"][arr["friction_type"] == 4] = 2

        # transform to a Pipes object
        return Pipes(**{name: arr[name] for name in arr.dtype.names})

    def get_pumps(self) -> Pumps:
        with self.get_session() as session:
            arr = (
                session.query(
                    models.Pumpstation.id,
                    models.Pumpstation.code,
                    models.Pumpstation.capacity,
                    models.Pumpstation.connection_node_start_id,
                    models.Pumpstation.connection_node_end_id,
                    models.Pumpstation.type_,
                    models.Pumpstation.start_level,
                    models.Pumpstation.lower_stop_level,
                    models.Pumpstation.upper_stop_level,
                    models.Pumpstation.display_name,
                    models.Pumpstation.zoom_category,
                )
                .order_by(models.Pumpstation.id)
                .as_structarray()
            )

        # Pump capicity is entered as L/s but we need m3/s.
        arr["capacity"] = arr["capacity"] / 1000

        # transform to a Pumps object
        return Pumps(**{name: arr[name] for name in arr.dtype.names})

    def get_weirs(self) -> Weirs:
        """Return Weirs"""
        with self.get_session() as session:
            arr = (
                session.query(
                    models.Weir.id,
                    models.Weir.code,
                    models.Weir.connection_node_start_id,
                    models.Weir.connection_node_end_id,
                    models.Weir.crest_level,
                    models.Weir.crest_type,
                    models.Weir.cross_section_definition_id,
                    models.Weir.discharge_coefficient_negative,
                    models.Weir.discharge_coefficient_positive,
                    models.Weir.friction_type,
                    models.Weir.friction_value,
                    models.Weir.display_name,
                    models.Weir.zoom_category,
                    models.Weir.sewerage,
                )
                .order_by(models.Weir.id)
                .as_structarray()
            )

        # map friction_type 4 to friction_type 2 to match crosssectionlocation enum
        arr["friction_type"][arr["friction_type"] == 4] = 2

        return Weirs(**{name: arr[name] for name in arr.dtype.names})


# Constructing a Transformer takes quite long, so we use caching here. The
# function is deterministic so this doesn't have any side effects.
@lru_cache(maxsize=8)
def _get_reproject_func(source_epsg: int, target_epsg: int) -> Callable:
    transformer = Transformer.from_crs(
        CRS.from_epsg(source_epsg), CRS.from_epsg(target_epsg), always_xy=True
    )

    def func(coords):
        if coords.shape[0] == 0:
            return coords
        x, y = transformer.transform(coords[:, 0], coords[:, 1])
        return np.array([x, y]).T

    return func
