# Module containing functions to read from the SQLite

from condenser import NumpyQuery
from contextlib import contextmanager
from functools import lru_cache
from pyproj import Transformer
from pyproj.crs import CRS
from sqlalchemy import cast
from sqlalchemy import inspect
from sqlalchemy import Integer
from sqlalchemy.orm import Session
from threedi_modelchecker.threedi_database import ThreediDatabase
from threedi_modelchecker.threedi_model import models
from threedi_modelchecker.threedi_model.custom_types import IntegerEnum
from threedigrid_builder.base import MakeGridSettings, MakeTablesSettings
from threedigrid_builder.grid import Channels, GridAttrs
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import CrossSectionDefinitions
from threedigrid_builder.grid import CrossSectionLocations
from threedigrid_builder.grid import Culverts
from threedigrid_builder.grid import GridRefinements
from threedigrid_builder.grid import Orifices
from threedigrid_builder.grid import Pipes
from threedigrid_builder.grid import Weirs
from typing import Callable, Tuple
from typing import ContextManager

import numpy as np
import pathlib
import pygeos


__all__ = ["SQLite"]


# hardcoded source projection
SOURCE_EPSG = 4326

# put some global defaults on datatypes
NumpyQuery.default_numpy_settings[Integer] = {"dtype": np.int32, "null": -9999}
NumpyQuery.default_numpy_settings[IntegerEnum] = {
    **NumpyQuery.default_numpy_settings[Integer],
    "sql_cast": lambda x: cast(x, Integer),
}


class SQLite:
    def __init__(self, path: pathlib.Path):
        path = str(path)
        sqlite_settings = {"db_path": path, "db_file": path}
        self.db = ThreediDatabase(
            connection_settings=sqlite_settings, db_type="spatialite"
        )
        self._epsg_code = None  # for reproject()

    @property
    def epsg_code(self) -> int:
        if self._epsg_code is None:
            self._epsg_code = self.get_settings()[0]
        return self._epsg_code

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
           dict with epsg_code, model_name, make_grid_settings, make_table_settings
        """
        with self.get_session() as session:
            global_ = _object_as_dict(session.query(models.GlobalSetting).order_by("id").first())
            groundwater = _object_as_dict(session.query(models.GroundWater).order_by("id").first())
            interflow = _object_as_dict(session.query(models.Interflow).order_by("id").first())
            infiltration = _object_as_dict(session.query(models.SimpleInfiltration).order_by("id").first())

        make_grid = MakeGridSettings.from_dict(global_)
        make_tables = MakeTablesSettings.from_dict(
            {**groundwater, **interflow, **infiltration, **global_}
        )
        return {
            "epsg_code": global_["epsg_code"],
            "model_name": global_["name"],
            "make_grid_settings": make_grid,
            "make_tables_settings": make_tables,
        }

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
                    models.Manhole.id.label("manhole_id"),
                    models.Manhole.calculation_type,
                    models.Manhole.bottom_level,
                    models.Manhole.drain_level,
                    models.Manhole.manhole_indicator,
                    models.Manhole.surface_level,
                    models.Manhole.shape.label("manhole_shape"),
                    models.Manhole.width.label("manhole_width"),
                )
                .join(models.ConnectionNode.manholes, isouter=True)
                .distinct(models.ConnectionNode.id)
                .order_by(models.ConnectionNode.id)
                .as_structarray()
            )

        arr["the_geom"] = self.reproject(arr["the_geom"])

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
                )
                .order_by(models.Culvert.id)
                .as_structarray()
            )

        arr["the_geom"] = self.reproject(arr["the_geom"])

        # map friction_type 4 to friction_type 2 to match crosssectionlocation enum
        arr["friction_type"][arr["friction_type"] == 4] = 2

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
                .order_by(models.GridRefinementArea.id)
                .as_structarray()
            )
            arr = np.concatenate((arr1, arr2))

        # reproject
        arr["the_geom"] = self.reproject(arr["the_geom"])
        arr["id"] = np.arange(len(arr["refinement_level"]))

        return GridRefinements(**{name: arr[name] for name in arr.dtype.names})

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
                )
                .order_by(models.Pipe.id)
                .as_structarray()
            )

        # map friction_type 4 to friction_type 2 to match crosssectionlocation enum
        arr["friction_type"][arr["friction_type"] == 4] = 2

        # transform to a Pipes object
        return Pipes(**{name: arr[name] for name in arr.dtype.names})

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
                )
                .order_by(models.Weir.id)
                .as_structarray()
            )

        # map friction_type 4 to friction_type 2 to match crosssectionlocation enum
        arr["friction_type"][arr["friction_type"] == 4] = 2

        return Weirs(**{name: arr[name] for name in arr.dtype.names})


def _object_as_dict(obj) -> dict:
    if obj is None:
        return {}
    # https://stackoverflow.com/questions/1958219/convert-sqlalchemy-row-object-to-python-dict
    return {c.key: getattr(obj, c.key) for c in inspect(obj).mapper.column_attrs}


# Constructing a Transformer takes quite long, so we use caching here. The
# function is deterministic so this doesn't have any side effects.
@lru_cache(maxsize=8)
def _get_reproject_func(source_epsg: int, target_epsg: int) -> Callable:
    transformer = Transformer.from_crs(
        CRS.from_epsg(source_epsg), CRS.from_epsg(target_epsg), always_xy=True
    )

    def func(coords):
        x, y = transformer.transform(coords[:, 0], coords[:, 1])
        return np.array([x, y]).T

    return func
