# Module containing functions to read from the SQLite

import pathlib
from contextlib import contextmanager
from enum import Enum
from typing import ContextManager, Dict, List, Union

import numpy as np
import shapely
from condenser import NumpyQuery
from sqlalchemy import case, cast, func, inspect, Integer, literal
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

MIN_SQLITE_VERSION = 229

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


def map_cross_section_definition(
    objects: List[Union[CrossSectionLocations, Pipes, Weirs, Orifices, Culverts]],
    definition_map: Dict[str, Dict[int, int]],
) -> None:
    """
    Set cross section definition ids for cross_section_locations,
    pipes, weirs, orifices and culverts to match the unique
    cross section locations.

    Args:
        objects: List of objects (CrossSectionLocations, Pipes, Weirs, Orifices, Culverts) to map the definition to
        definition_map: Mapping of object names to their definition IDs
    """
    # Map shapes to key in definition_map
    object_map = {
        CrossSectionLocations: "cross_section_location",
        Pipes: "pipe",
        Weirs: "weir",
        Orifices: "orifice",
        Culverts: "culvert",
    }
    for object in objects:
        object_name = object_map.get(type(object))
        if object_name is None:
            raise ValueError(f"Object of type {type(object)} cannot be mapped")
        mapping = definition_map.get(object_name)
        if mapping is None:
            continue
        idx = object.id_to_index(list(mapping.keys()))
        # set correct cross section definition id's for mapped object
        if isinstance(object, CrossSectionLocations):
            object.definition_id[idx] = np.array(list(mapping.values()), dtype=int)
        else:
            object.cross_section_definition_id[idx] = np.array(
                list(mapping.values()), dtype=int
            )


class SQLite:
    def __init__(self, path: pathlib.Path, upgrade=False, convert_to_geopackage=False):
        if not path.exists():
            raise FileNotFoundError(f"File not found: {path}")
        self.db = ThreediDatabase(path)
        self._epsg_code = None
        version = self.get_version()
        if version < MIN_SQLITE_VERSION:
            if upgrade:
                self.upgrade(convert_to_geopackage=convert_to_geopackage)
            else:
                raise SchematisationError(f"Too old sqlite version {version}.")

    @property
    def epsg_code(self) -> int:
        if self._epsg_code is None:
            self._epsg_code = self.db.schema.epsg_code
        return self._epsg_code

    @epsg_code.setter
    def epsg_code(self, value: int):
        self._epsg_code = value

    def get_version(self) -> int:
        # check version
        schema = ModelSchema(self.db)
        return schema.get_version()

    def upgrade(self, convert_to_geopackage=False):
        schema = ModelSchema(self.db)
        schema.upgrade(backup=False)

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
                    session.query(models.VegetationDrag2D).one()
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
        if interception:
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
            _set_initialization_type(
                groundwater,
                "equilibrium_infiltration_rate",
                type_field="equilibrium_infiltration_rate_aggregation",
            )
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
            "model_name": model_settings["name"],
            "grid_settings": grid_settings,
            "tables_settings": tables_settings,
        }

    def get_surfaces(self) -> Surfaces:
        with self.get_session() as session:
            arr = (
                session.query(
                    models.Surface.id.label("surface_id"),
                    models.Surface.code,
                    models.Surface.display_name,
                    models.Surface.area,
                    models.Surface.geom.label("the_geom"),
                    models.SurfaceParameters.outflow_delay,
                    models.SurfaceParameters.surface_layer_thickness,
                    models.SurfaceParameters.infiltration,
                    models.SurfaceParameters.max_infiltration_capacity,
                    models.SurfaceParameters.min_infiltration_capacity,
                    models.SurfaceParameters.infiltration_decay_constant,
                    models.SurfaceParameters.infiltration_recovery_constant,
                    models.SurfaceMap.connection_node_id,
                    models.SurfaceMap.percentage,
                    models.ConnectionNode.geom.label("connection_node_the_geom"),
                )
                .select_from(models.Surface)
                .join(
                    models.SurfaceParameters,
                    models.Surface.surface_parameters_id == models.SurfaceParameters.id,
                )
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
                    models.BoundaryCondition1D.type,
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
                    models.BoundaryConditions2D.type,
                    models.BoundaryConditions2D.geom.label("the_geom"),
                )
                .order_by(models.BoundaryConditions2D.id)
                .as_structarray()
            )
        # transform to a BoundaryConditions1D object
        return BoundaryConditions2D(**{name: arr[name] for name in arr.dtype.names})

    def get_channels(self) -> Channels:
        """Return Channels"""
        cols = [
            models.Channel.geom.label("the_geom"),
            models.Channel.calculation_point_distance.label("dist_calc_points"),
            models.Channel.id,
            models.Channel.code,
            models.Channel.connection_node_id_start.label("connection_node_start_id"),
            models.Channel.connection_node_id_end.label("connection_node_end_id"),
            models.Channel.exchange_type.label("calculation_type"),
            models.Channel.display_name,
            models.Channel.exchange_thickness,
            models.Channel.hydraulic_conductivity_out,
            models.Channel.hydraulic_conductivity_in,
        ]

        with self.get_session() as session:
            arr = session.query(*cols).order_by(models.Channel.id).as_structarray()

        # map "old" calculation types (100, 101, 102, 105) to (0, 1, 2, 5)
        arr["calculation_type"][arr["calculation_type"] >= 100] -= 100
        arr["hydraulic_conductivity_out"] /= DAY_IN_SECONDS
        arr["hydraulic_conductivity_in"] /= DAY_IN_SECONDS

        # transform to a Channels object
        return Channels(**{name: arr[name] for name in arr.dtype.names})

    def get_connection_nodes(self) -> ConnectionNodes:
        """Return ConnectionNodes"""
        cols = [
            models.ConnectionNode.geom.label("the_geom"),
            models.ConnectionNode.id,
            models.ConnectionNode.code,
            models.ConnectionNode.storage_area,
            models.ConnectionNode.initial_water_level.label("initial_waterlevel"),
            models.ConnectionNode.exchange_type.label("calculation_type"),
            models.ConnectionNode.bottom_level,
            models.ConnectionNode.exchange_level.label("drain_level"),
            models.ConnectionNode.visualisation.label("manhole_indicator"),
            models.ConnectionNode.manhole_surface_level.label("surface_level"),
            models.ConnectionNode.display_name,
            models.ConnectionNode.exchange_thickness,
            models.ConnectionNode.hydraulic_conductivity_out,
            models.ConnectionNode.hydraulic_conductivity_in,
        ]

        with self.get_session() as session:
            arr = (
                session.query(*cols).order_by(models.ConnectionNode.id).as_structarray()
            )

        # replace -9999.0 with NaN in initial_water_level
        arr["initial_waterlevel"][arr["initial_waterlevel"] == -9999.0] = np.nan
        arr["hydraulic_conductivity_out"] /= DAY_IN_SECONDS
        arr["hydraulic_conductivity_in"] /= DAY_IN_SECONDS

        return ConnectionNodes(**{name: arr[name] for name in arr.dtype.names})

    def get_cross_section_definition_for_table(self, table) -> np.ndarray:
        with self.get_session() as session:
            cols = [
                literal(table.__tablename__).label("origin_table"),
                table.id.label("origin_id"),
                table.cross_section_shape.label("shape"),
                table.cross_section_width.label("width"),
                table.cross_section_height.label("height"),
                table.cross_section_table,
            ]
            if table == models.CrossSectionLocation:
                cols += [
                    table.cross_section_friction_values.label("friction_values"),
                    table.cross_section_vegetation_table,
                    table.vegetation_stem_density,
                    table.vegetation_stem_diameter,
                    table.vegetation_height,
                    table.vegetation_drag_coefficient,
                ]
            arr = session.query(*cols).select_from(table).as_structarray()
            # map shape 10 to 1 (circle) to match CrossSectionShape enum
            arr["shape"][arr["shape"] == 10] = 1
            # map shape 11 to 5 (tabulated rectangle) to match CrossSectionShape enum
            arr["shape"][arr["shape"] == 11] = 5
            # map shape 12 to 6 (tabulated trapezium) to match CrossSectionShape enum
            arr["shape"][arr["shape"] == 12] = 6
        return arr

    def get_cross_section_definitions(self) -> CrossSectionDefinitions:
        """Return CrossSectionDefinitions"""
        attr_dict = {
            "id": np.empty(0, dtype="i4"),
            "origin_table": np.empty(0, dtype="O"),
            "origin_id": np.empty(0, dtype="i4"),
            "shape": np.empty(0, dtype="O"),
            "width": np.empty(0, dtype="f8"),
            "height": np.empty(0, dtype="f8"),
            "cross_section_table": np.empty(0, dtype="O"),
            "friction_values": np.empty(0, dtype="O"),
            "cross_section_vegetation_table": np.empty(0, dtype="O"),
        }
        for table in [
            models.CrossSectionLocation,
            models.Culvert,
            models.Orifice,
            models.Pipe,
            models.Weir,
        ]:
            arr = self.get_cross_section_definition_for_table(table)
            if len(arr) == 0:
                continue
            for name in attr_dict.keys():
                if name == "id":
                    data = len(attr_dict["id"]) + np.arange(len(arr))
                elif name in arr.dtype.names:
                    data = arr[name]
                else:
                    data = np.empty(len(arr), attr_dict[name].dtype)
                attr_dict[name] = np.concatenate((attr_dict[name], data))

        return CrossSectionDefinitions(**attr_dict)

    def get_cross_section_locations(self) -> CrossSectionLocations:
        """Return CrossSectionLocations"""
        with self.get_session() as session:
            arr = (
                session.query(
                    models.CrossSectionLocation.id,
                    models.CrossSectionLocation.code,
                    models.CrossSectionLocation.geom.label("the_geom"),
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

        # transform to a CrossSectionLocations object
        return CrossSectionLocations(**{name: arr[name] for name in arr.dtype.names})

    def get_culverts(self) -> Culverts:
        """Return Culverts"""
        cols = [
            models.Culvert.id,
            models.Culvert.code,
            models.Culvert.geom.label("the_geom"),
            models.Culvert.calculation_point_distance.label("dist_calc_points"),
            models.Culvert.connection_node_id_start.label("connection_node_start_id"),
            models.Culvert.connection_node_id_end.label("connection_node_end_id"),
            models.Culvert.exchange_type.label("calculation_type"),
            models.Culvert.invert_level_start.label("invert_level_start_point"),
            models.Culvert.invert_level_end.label("invert_level_end_point"),
            models.Culvert.discharge_coefficient_negative,
            models.Culvert.discharge_coefficient_positive,
            models.Culvert.display_name,
            case(
                {
                    models.Culvert.friction_value.isnot(None)
                    & models.Culvert.friction_type.isnot(
                        None
                    ): models.Culvert.friction_value
                },
                else_=models.Material.friction_coefficient,
            ).label("friction_value"),
            case(
                {
                    models.Culvert.friction_value.isnot(None)
                    & models.Culvert.friction_type.isnot(
                        None
                    ): models.Culvert.friction_type
                },
                else_=models.Material.friction_type,
            ).label("friction_type"),
        ]

        with self.get_session() as session:
            arr = (
                session.query(*cols)
                .outerjoin(
                    models.Material, models.Culvert.material_id == models.Material.id
                )
                .order_by(models.Culvert.id)
                .as_structarray()
            )

        # map friction_type 4 to friction_type 2 to match crosssectionlocation enum
        arr["friction_type"][arr["friction_type"] == 4] = 2

        # When no calculation type is provides we default to isolated
        arr["calculation_type"][
            arr["calculation_type"] == -9999
        ] = LineType.LINE_1D_ISOLATED
        # map "old" calculation types (100, 101, 102, 105) to (0, 1, 2, 5)
        arr["calculation_type"][arr["calculation_type"] >= 100] -= 100

        # transform to a Culverts object
        return Culverts(**{name: arr[name] for name in arr.dtype.names})

    def get_exchange_lines(self) -> ExchangeLines:
        with self.get_session() as session:
            arr = (
                session.query(
                    models.ExchangeLine.id,
                    models.ExchangeLine.channel_id,
                    models.ExchangeLine.geom.label("the_geom"),
                    models.ExchangeLine.exchange_level,
                )
                .order_by(models.ExchangeLine.id)
                .as_structarray()
            )

        # transform to a Channels object
        return ExchangeLines(**{name: arr[name] for name in arr.dtype.names})

    def get_grid_refinements(self) -> GridRefinements:
        """Return Gridrefinement and GridRefinementArea concatenated into one array."""
        with self.get_session() as session:
            arr1 = (
                session.query(
                    models.GridRefinementLine.geom.label("the_geom"),
                    models.GridRefinementLine.id,
                    models.GridRefinementLine.code,
                    models.GridRefinementLine.display_name,
                    models.GridRefinementLine.grid_level.label("refinement_level"),
                )
                .filter(
                    models.GridRefinementLine.geom.isnot(None),
                    models.GridRefinementLine.grid_level.isnot(None),
                )
                .order_by(models.GridRefinementLine.id)
                .as_structarray()
            )
            arr2 = (
                session.query(
                    models.GridRefinementArea.geom.label("the_geom"),
                    models.GridRefinementArea.id,
                    models.GridRefinementArea.code,
                    models.GridRefinementArea.display_name,
                    models.GridRefinementArea.grid_level.label("refinement_level"),
                )
                .filter(
                    models.GridRefinementArea.geom.isnot(None),
                    models.GridRefinementArea.grid_level.isnot(None),
                )
                .order_by(models.GridRefinementArea.id)
                .as_structarray()
            )
            arr = np.concatenate((arr1, arr2))

        arr["id"] = np.arange(len(arr["refinement_level"]))
        return GridRefinements(**{name: arr[name] for name in arr.dtype.names})

    def get_dem_average_areas(self) -> DemAverageAreas:
        """Return DemAverageAreas"""
        with self.get_session() as session:
            arr = (
                session.query(
                    models.DemAverageArea.id,
                    models.DemAverageArea.geom.label("the_geom"),
                )
                .order_by(models.DemAverageArea.id)
                .as_structarray()
            )

        return DemAverageAreas(**{name: arr[name] for name in arr.dtype.names})

    def get_obstacles(self) -> Obstacles:
        with self.get_session() as session:
            arr = (
                session.query(
                    models.Obstacle.geom.label("the_geom"),
                    models.Obstacle.id,
                    models.Obstacle.crest_level,
                    models.Obstacle.affects_2d,
                    models.Obstacle.affects_1d2d_closed,
                    models.Obstacle.affects_1d2d_open_water,
                )
                .order_by(models.Obstacle.id)
                .as_structarray()
            )
        return Obstacles(**{name: arr[name] for name in arr.dtype.names})

    def get_orifices(self) -> Orifices:
        """Return Orifices"""
        cols = [
            models.Orifice.id,
            models.Orifice.code,
            models.Orifice.connection_node_id_start.label("connection_node_start_id"),
            models.Orifice.connection_node_id_end.label("connection_node_end_id"),
            models.Orifice.crest_level,
            models.Orifice.crest_type,
            models.Orifice.discharge_coefficient_negative,
            models.Orifice.discharge_coefficient_positive,
            models.Orifice.display_name,
            models.Orifice.sewerage,
            models.Orifice.geom.label("the_geom"),
            case(
                {
                    models.Orifice.friction_value.isnot(None)
                    & models.Orifice.friction_type.isnot(
                        None
                    ): models.Orifice.friction_value
                },
                else_=models.Material.friction_coefficient,
            ).label("friction_value"),
            case(
                {
                    models.Orifice.friction_value.isnot(None)
                    & models.Orifice.friction_type.isnot(
                        None
                    ): models.Orifice.friction_type
                },
                else_=models.Material.friction_type,
            ).label("friction_type"),
        ]
        with self.get_session() as session:
            arr = (
                session.query(*cols)
                .outerjoin(
                    models.Material, models.Orifice.material_id == models.Material.id
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
            models.Pipe.exchange_type.label("calculation_type"),
            models.Pipe.invert_level_start.label("invert_level_start_point"),
            models.Pipe.invert_level_end.label("invert_level_end_point"),
            models.Pipe.calculation_point_distance.label("dist_calc_points"),
            models.Pipe.connection_node_id_start.label("connection_node_start_id"),
            models.Pipe.connection_node_id_end.label("connection_node_end_id"),
            models.Pipe.display_name,
            models.Pipe.exchange_thickness,
            models.Pipe.hydraulic_conductivity_out,
            models.Pipe.hydraulic_conductivity_in,
            models.Pipe.material_id.label("material"),
            models.Pipe.geom.label("the_geom"),
            case(
                {
                    models.Pipe.friction_value.isnot(None)
                    & models.Pipe.friction_type.isnot(None): models.Pipe.friction_value
                },
                else_=models.Material.friction_coefficient,
            ).label("friction_value"),
            case(
                {
                    models.Pipe.friction_value.isnot(None)
                    & models.Pipe.friction_type.isnot(None): models.Pipe.friction_type
                },
                else_=models.Material.friction_type,
            ).label("friction_type"),
        ]
        with self.get_session() as session:
            arr = (
                session.query(*cols)
                .outerjoin(
                    models.Material, models.Pipe.material_id == models.Material.id
                )
                .order_by(models.Pipe.id)
                .as_structarray()
            )

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
                    models.Pump.id,
                    models.Pump.code,
                    models.Pump.capacity,
                    models.Pump.connection_node_id.label("connection_node_start_id"),
                    models.PumpMap.connection_node_id_end.label(
                        "connection_node_end_id"
                    ),
                    models.Pump.type_,
                    models.Pump.start_level,
                    models.Pump.lower_stop_level,
                    models.Pump.upper_stop_level,
                    models.Pump.display_name,
                )
                .outerjoin(models.PumpMap, models.Pump.id == models.PumpMap.pump_id)
                .order_by(models.Pump.id)
                .as_structarray()
            )

        # Pump capacity is entered as L/s but we need m3/s.
        arr["capacity"] = arr["capacity"] / 1000

        # transform to a Pumps object
        return Pumps(**{name: arr[name] for name in arr.dtype.names})

    def get_weirs(self) -> Weirs:
        """Return Weirs"""
        cols = [
            models.Weir.id,
            models.Weir.code,
            models.Weir.connection_node_id_start.label("connection_node_start_id"),
            models.Weir.connection_node_id_end.label("connection_node_end_id"),
            models.Weir.crest_level,
            models.Weir.crest_type,
            models.Weir.discharge_coefficient_negative,
            models.Weir.discharge_coefficient_positive,
            models.Weir.display_name,
            models.Weir.sewerage,
            models.Weir.geom.label("the_geom"),
            case(
                {
                    models.Weir.friction_value.isnot(None)
                    & models.Weir.friction_type.isnot(None): models.Weir.friction_value
                },
                else_=models.Material.friction_coefficient,
            ).label("friction_value"),
            case(
                {
                    models.Weir.friction_value.isnot(None)
                    & models.Weir.friction_type.isnot(None): models.Weir.friction_type
                },
                else_=models.Material.friction_type,
            ).label("friction_type"),
        ]
        with self.get_session() as session:
            arr = (
                session.query(*cols)
                .outerjoin(
                    models.Material, models.Weir.material_id == models.Material.id
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
                    models.Windshielding1D.id,
                    models.Windshielding1D.channel_id,
                    models.Windshielding1D.north,
                    models.Windshielding1D.northeast,
                    models.Windshielding1D.east,
                    models.Windshielding1D.southeast,
                    models.Windshielding1D.south,
                    models.Windshielding1D.southwest,
                    models.Windshielding1D.west,
                    models.Windshielding1D.northwest,
                )
                .order_by(models.Windshielding1D.id)
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
            models.PotentialBreach.geom.label("the_geom"),
            models.PotentialBreach.channel_id,
        ]

        if self.get_version() >= 212:
            cols += [
                models.PotentialBreach.initial_exchange_level.label("exchange_level"),
                models.PotentialBreach.final_exchange_level.label(
                    "maximum_breach_depth"
                ),
                models.PotentialBreach.levee_material,
            ]

        with self.get_session() as session:
            arr = (
                session.query(*cols)
                .order_by(models.PotentialBreach.id)
                .as_structarray()
            )

        # derive maximum_breach_depth from initial and final exchange level
        # and overwrite final_exchange_level because adding a field is more work
        arr["maximum_breach_depth"] = (
            arr["exchange_level"] - arr["maximum_breach_depth"]
        )

        return PotentialBreaches(**{name: arr[name] for name in arr.dtype.names})
