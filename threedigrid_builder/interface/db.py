# Module containing functions to read from the SQLite

import pathlib
from contextlib import contextmanager
from enum import Enum
from functools import lru_cache
from typing import Callable, ContextManager

import numpy as np
import shapely
from condenser import NumpyQuery
from pyproj import Transformer
from pyproj.crs import CRS
from sqlalchemy import cast, func, inspect, Integer
from sqlalchemy.orm import Session
from threedi_schema import custom_types, models, ModelSchema, ThreediDatabase

from threedigrid_builder.base import GridSettings, Pumps, TablesSettings
from threedigrid_builder.constants import InitializationType, LineType
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import (
    BoundaryConditions1D,
    BoundaryConditions2D,
    Channels,
    ConnectionNodes,
    CrossSectionDefinitions,
    CrossSectionLocations,
    Culverts,
    DemAverageAreas,
    ExchangeLines,
    GridRefinements,
    Obstacles,
    Orifices,
    Pipes,
    PotentialBreaches,
    Surfaces,
    Weirs,
    Windshieldings,
)

__all__ = ["SQLite"]

# hardcoded source projection
SOURCE_EPSG = 4326

MIN_SQLITE_VERSION = 223

DAY_IN_SECONDS = 24.0 * 3600.0

# put some global defaults on datatypes
NumpyQuery.default_numpy_settings[Integer] = {"dtype": np.int32, "null": -9999}
NumpyQuery.default_numpy_settings[custom_types.IntegerEnum] = {
    **NumpyQuery.default_numpy_settings[Integer],
    "sql_cast": lambda x: cast(x, Integer),
}
NumpyQuery.default_numpy_settings[custom_types.Geometry] = {
    "dtype": np.dtype("O"),
    "sql_cast": func.ST_AsBinary,
    "numpy_cast": shapely.from_wkb,
}


def _object_as_dict(obj) -> dict:
    """Convert SQLAlchemy object to dict, casting Enums. Optional to prefix keys."""
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
    def __init__(self, path: pathlib.Path, upgrade=False, convert_to_geopackage=False):
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")
        self.db = ThreediDatabase(path)
        self._epsg_code = None  # for reproject()

        version = self.get_version()
        if version < MIN_SQLITE_VERSION:
            if upgrade:
                self.upgrade(convert_to_geopackage=convert_to_geopackage)
            else:
                raise SchematisationError(f"Too old sqlite version {version}.")

    def get_version(self) -> int:
        # check version
        schema = ModelSchema(self.db)
        return schema.get_version()

    def upgrade(self, convert_to_geopackage=False):
        schema = ModelSchema(self.db)
        schema.upgrade(
            backup=False, set_views=False, convert_to_geopackage=convert_to_geopackage
        )

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
            model_settings = session.query(models.ModelSettings).order_by("id").first()
            if (
                model_settings.use_groundwater_flow
                or model_settings.use_groundwater_storage
            ):
                groundwater = _object_as_dict(session.query(models.GroundWater).one())
            else:
                groundwater = {}
            if model_settings.use_interflow:
                interflow = _object_as_dict(session.query(models.Interflow).one())
            else:
                interflow = {}
            if model_settings.use_simple_infiltration:
                infiltration = _object_as_dict(
                    session.query(models.SimpleInfiltration).one()
                )
            else:
                infiltration = {}
            if model_settings.use_vegetation_drag_2d:
                vegetation_drag = _object_as_dict(
                    session.query(models.VegetationDrag).one()
                )
            else:
                vegetation_drag = {}
            if model_settings.use_interception:
                interception = _object_as_dict(session.query(models.Interception).one())
            else:
                interception = {}
            model_settings = _object_as_dict(model_settings)

        # record if there is a DEM file to be expected
        # Note: use_2d_flow only determines whether there are flow lines
        model_settings["use_2d"] = bool(model_settings["dem_file"])

        # set/adapt initialization types to include information about file presence
        NO_AGG = InitializationType.NO_AGG
        AVERAGE = InitializationType.AVERAGE
        _set_initialization_type(
            model_settings,
            "friction_coefficient",
            default=AVERAGE if model_settings["friction_averaging"] else NO_AGG,
        )
        _set_initialization_type(
            interception,
            "interception",
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
            _set_initialization_type(
                infiltration, "max_infiltration_volume", default=NO_AGG
            )
        if groundwater:
            # default is what the user supplied (MIN/MAX/AVERAGE)

            _set_initialization_type(
                groundwater,
                "groundwater_impervious_layer_level",
                type_field="groundwater_impervious_layer_level_aggregation",
            )
            _set_initialization_type(
                groundwater,
                "phreatic_storage_capacity",
                type_field="phreatic_storage_capacity_aggregation",
            )
            _set_initialization_type(groundwater, "equilibrium_infiltration_rate")
            _set_initialization_type(
                groundwater,
                "initial_infiltration_rate",
                type_field="initial_infiltration_rate_aggregation",
            )
            _set_initialization_type(
                groundwater,
                "infiltration_decay_period",
                type_field="infiltration_decay_period_aggregation",
            )
            _set_initialization_type(
                groundwater,
                "groundwater_hydraulic_conductivity",
                type_field="groundwater_hydraulic_conductivity_aggregation",
            )

        if vegetation_drag:
            _set_initialization_type(
                vegetation_drag, "vegetation_height", default=NO_AGG
            )
            _set_initialization_type(
                vegetation_drag, "vegetation_stem_count", default=NO_AGG
            )
            _set_initialization_type(
                vegetation_drag, "vegetation_stem_diameter", default=NO_AGG
            )
            _set_initialization_type(
                vegetation_drag, "vegetation_drag_coefficient", default=NO_AGG
            )
        # Copy Simulation Template Settings to model_settings dict
        template_settings = _object_as_dict(
            session.query(models.SimulationTemplateSettings).order_by("id").first()
        )
        model_settings["name"] = template_settings["name"]
        model_settings["use_0d_inflow"] = template_settings["use_0d_inflow"]

        grid_settings = GridSettings.from_dict(model_settings)
        tables_settings = TablesSettings.from_dict(
            {
                **groundwater,
                **interflow,
                **infiltration,
                **vegetation_drag,
                **model_settings,
                **interception,
            }
        )
        return {
            "epsg_code": model_settings["epsg_code"],
            "model_name": model_settings["name"],
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
          shapely+pyproj is approx 2x faster than spatialite

        Args:
          geometries (ndarray of shapely.Geometry): geometries in EPSG 4326
        """
        target_epsg = self.epsg_code
        func = _get_reproject_func(SOURCE_EPSG, target_epsg)
        return shapely.transform(geometries, func)

    def get_surfaces(self) -> Surfaces:
        with self.get_session() as session:
            arr = (
                session.query(
                    models.Surface.id.label("surface_id"),
                    models.Surface.code,
                    models.Surface.display_name,
                    models.Surface.area,
                    models.Surface.geom,
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
        arr["geom"] = self.reproject(arr["geom"])
        arr["connection_node_the_geom"] = self.reproject(
            arr["connection_node_the_geom"]
        )

        return Surfaces(
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

    def get_channels(self) -> Channels:
        """Return Channels"""
        cols = [
            models.Channel.the_geom,
            models.Channel.dist_calc_points,
            models.Channel.id,
            models.Channel.code,
            models.Channel.connection_node_start_id,
            models.Channel.connection_node_end_id,
            models.Channel.calculation_type,
            models.Channel.display_name,
            models.Channel.zoom_category,
            models.Channel.exchange_thickness,
            models.Channel.hydraulic_conductivity_out,
            models.Channel.hydraulic_conductivity_in,
        ]

        with self.get_session() as session:
            arr = session.query(*cols).order_by(models.Channel.id).as_structarray()

        arr["the_geom"] = self.reproject(arr["the_geom"])
        # map "old" calculation types (100, 101, 102, 105) to (0, 1, 2, 5)
        arr["calculation_type"][arr["calculation_type"] >= 100] -= 100
        arr["hydraulic_conductivity_out"] /= DAY_IN_SECONDS
        arr["hydraulic_conductivity_in"] /= DAY_IN_SECONDS

        # transform to a Channels object
        return Channels(**{name: arr[name] for name in arr.dtype.names})

    def get_connection_nodes(self) -> ConnectionNodes:
        """Return ConnectionNodes (which are enriched using the manhole table)"""
        cols = [
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
            models.Manhole.exchange_thickness,
            models.Manhole.hydraulic_conductivity_out,
            models.Manhole.hydraulic_conductivity_in,
        ]

        with self.get_session() as session:
            arr = (
                session.query(*cols)
                .join(models.ConnectionNode.manholes, isouter=True)
                .order_by(models.ConnectionNode.id)
                .as_structarray()
            )

        arr["the_geom"] = self.reproject(arr["the_geom"])

        # replace -9999.0 with NaN in initial_waterlevel
        arr["initial_waterlevel"][arr["initial_waterlevel"] == -9999.0] = np.nan
        arr["hydraulic_conductivity_out"] /= DAY_IN_SECONDS
        arr["hydraulic_conductivity_in"] /= DAY_IN_SECONDS

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
                    models.CrossSectionDefinition.friction_values,
                    models.CrossSectionDefinition.vegetation_stem_densities,
                    models.CrossSectionDefinition.vegetation_stem_diameters,
                    models.CrossSectionDefinition.vegetation_heights,
                    models.CrossSectionDefinition.vegetation_drag_coefficients,
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
                    models.CrossSectionLocation.vegetation_stem_density,
                    models.CrossSectionLocation.vegetation_stem_diameter,
                    models.CrossSectionLocation.vegetation_height,
                    models.CrossSectionLocation.vegetation_drag_coefficient,
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

    def get_exchange_lines(self) -> ExchangeLines:
        with self.get_session() as session:
            arr = (
                session.query(
                    models.ExchangeLine.id,
                    models.ExchangeLine.channel_id,
                    models.ExchangeLine.the_geom,
                    models.ExchangeLine.exchange_level,
                )
                .order_by(models.ExchangeLine.id)
                .as_structarray()
            )

        arr["the_geom"] = self.reproject(arr["the_geom"])

        # transform to a Channels object
        return ExchangeLines(**{name: arr[name] for name in arr.dtype.names})

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
        cols = [
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
            models.Pipe.exchange_thickness,
            models.Pipe.hydraulic_conductivity_out,
            models.Pipe.hydraulic_conductivity_in,
        ]

        with self.get_session() as session:
            arr = session.query(*cols).order_by(models.Pipe.id).as_structarray()

        # map friction_type 4 to friction_type 2 to match crosssectionlocation enum
        arr["friction_type"][arr["friction_type"] == 4] = 2
        arr["hydraulic_conductivity_out"] /= DAY_IN_SECONDS
        arr["hydraulic_conductivity_in"] /= DAY_IN_SECONDS

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

    def get_windshieldings(self) -> Windshieldings:
        with self.get_session() as session:
            arr = (
                session.query(
                    models.Windshielding.id,
                    models.Windshielding.channel_id,
                    models.Windshielding.north,
                    models.Windshielding.northeast,
                    models.Windshielding.east,
                    models.Windshielding.southeast,
                    models.Windshielding.south,
                    models.Windshielding.southwest,
                    models.Windshielding.west,
                    models.Windshielding.northwest,
                )
                .order_by(models.Windshielding.id)
                .as_structarray()
            )

        return Windshieldings(**{name: arr[name] for name in arr.dtype.names})

    def get_potential_breaches(self) -> PotentialBreaches:
        # potential breaches are available from schema version 211
        if self.get_version() < 211:
            return PotentialBreaches(id=[])

        cols = [
            models.PotentialBreach.id,
            models.PotentialBreach.code,
            models.PotentialBreach.display_name,
            models.PotentialBreach.the_geom,
            models.PotentialBreach.channel_id,
        ]

        if self.get_version() >= 212:
            cols += [
                models.PotentialBreach.exchange_level,
                models.PotentialBreach.maximum_breach_depth,
                models.PotentialBreach.levee_material,
            ]

        with self.get_session() as session:
            arr = (
                session.query(*cols)
                .order_by(models.PotentialBreach.id)
                .as_structarray()
            )

        # reproject
        arr["the_geom"] = self.reproject(arr["the_geom"])

        return PotentialBreaches(**{name: arr[name] for name in arr.dtype.names})


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
