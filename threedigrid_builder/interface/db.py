# Module containing functions to read from the SQLite

from condenser import NumpyQuery
from contextlib import contextmanager
from functools import lru_cache
from pyproj import Transformer
from pyproj.crs import CRS
from sqlalchemy import cast
from sqlalchemy import inspect
from sqlalchemy import Integer
from threedi_modelchecker.threedi_database import ThreediDatabase
from threedi_modelchecker.threedi_model import models
from threedi_modelchecker.threedi_model.custom_types import IntegerEnum
from threedigrid_builder.grid import Channels
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import CrossSectionDefinitions
from threedigrid_builder.grid import CrossSectionLocations
from threedigrid_builder.grid import GridRefinements

import numpy as np
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
    def __init__(self, path):
        path = str(path)
        sqlite_settings = {"db_path": path, "db_file": path}
        self.db = ThreediDatabase(
            connection_settings=sqlite_settings, db_type="spatialite"
        )
        self._global_settings = None

    @contextmanager
    def get_session(self):
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

    @property
    def global_settings(self):
        """Return the global settings dictionary from the SQLite at path.

        The global settings are cached on self.
        """
        if self._global_settings is None:
            with self.get_session() as session:
                settings = session.query(models.GlobalSetting).order_by("id").first()
            self._global_settings = _object_as_dict(settings)
        return self._global_settings

    def reproject(self, geometries):
        """Reproject geometries from 4326 to the EPSG in the settings.

        Notes:
          pygeos+pyproj is approx 2x faster than spatialite

        Args:
          geometries (ndarray of pygeos.Geometry): geometries in EPSG 4326
        """
        target_epsg = self.global_settings["epsg_code"]
        func = _get_reproject_func(SOURCE_EPSG, target_epsg)
        return pygeos.apply(geometries, func)

    def get_channels(self):
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
        arr["calculation_type"] -= 100  # maps (100, 101, 102, 105) to (0, 1, 2, 5)

        # transform to a Channels object
        return Channels(**{name: arr[name] for name in arr.dtype.names})

    def get_connection_nodes(self):
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

    def get_cross_section_definitions(self):
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

    def get_cross_section_locations(self):
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

    def get_grid_refinements(self):
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


def _object_as_dict(obj):
    # https://stackoverflow.com/questions/1958219/convert-sqlalchemy-row-object-to-python-dict
    return {c.key: getattr(obj, c.key) for c in inspect(obj).mapper.column_attrs}


# Constructing a Transformer takes quite long, so we use caching here. The
# function is deterministic so this doesn't have any side effects.
@lru_cache(maxsize=8)
def _get_reproject_func(source_epsg, target_epsg):
    transformer = Transformer.from_crs(
        CRS.from_epsg(source_epsg), CRS.from_epsg(target_epsg), always_xy=True
    )

    def func(coords):
        x, y = transformer.transform(coords[:, 0], coords[:, 1])
        return np.array([x, y]).T

    return func
