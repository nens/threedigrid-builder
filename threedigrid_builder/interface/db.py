# Module containing functions to read from the SQLite

from contextlib import contextmanager
from functools import lru_cache
from geoalchemy2.functions import ST_AsBinary
from pyproj import Transformer
from pyproj.crs import CRS
from sqlalchemy import cast
from sqlalchemy import inspect
from sqlalchemy import Integer
from threedi_modelchecker.threedi_database import ThreediDatabase
from threedi_modelchecker.threedi_model import models
from threedigrid_builder.grid1d import Channels
from threedigrid_builder.grid1d import ConnectionNodes

import numpy as np
import pygeos


SOURCE_EPSG = 4326

__all__ = ["SQLite"]


class SQLite:
    def __init__(self, path):
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
        session = self.db.get_session()
        yield session
        session.close()

    @property
    def global_settings(self):
        """Return the global settings dictionary from the SQLite at path
        """
        if self._global_settings is None:
            with self.get_session() as session:
                settings = session.query(models.GlobalSetting).one()
            self._global_settings = _object_as_dict(settings)
        return self._global_settings

    def get_channels(self):
        """Return Channels

        Most stuff in this function will eventually be implemented in the package
        "condenser".
        """
        with self.get_session() as session:
            data = session.query(
                ST_AsBinary(models.Channel.the_geom),
                models.Channel.dist_calc_points,
                models.Channel.id,
                models.Channel.code,
                cast(models.Channel.connection_node_start_id, Integer),
                cast(models.Channel.connection_node_end_id, Integer),
                cast(models.Channel.calculation_type, Integer),
            ).all()

        # transform tuples to a numpy structured array
        arr = np.array(
            data,
            dtype=[
                ("the_geom", "O"),
                ("dist_calc_points", "f8"),
                ("id", "i8"),
                ("code", "O"),
                ("connection_node_start_id", "i8"),
                ("connection_node_end_id", "i8"),
                ("calculation_type", "i8"),
            ],
        )
        # transform to pygeos.Geometry
        arr["the_geom"] = pygeos.from_wkb(arr["the_geom"])

        # reproject
        target_epsg = self.global_settings["epsg_code"]
        arr["the_geom"] = pygeos.apply(
            arr["the_geom"], _get_reproject_func(SOURCE_EPSG, target_epsg)
        )

        # transform to a dict of 1D ndarrays
        return Channels(**{name: arr[name] for name in arr.dtype.names})

    def get_connection_nodes(self):
        """Return ConnectionNodes
        Most stuff in this function will eventually be implemented in the package
        "condenser".
        """
        with self.get_session() as session:
            data = session.query(
                ST_AsBinary(models.ConnectionNode.the_geom),
                models.ConnectionNode.id,
                models.ConnectionNode.code,
                models.ConnectionNode.storage_area,
            ).all()

        # transform tuples to a numpy structured array
        arr = np.array(
            data,
            dtype=[
                ("the_geom", "O"),
                ("id", "i8"),
                ("code", "O"),
                ("storage_area", "f8"),
            ],
        )
        # transform to pygeos.Geometry
        arr["the_geom"] = pygeos.from_wkb(arr["the_geom"])

        # reproject
        target_epsg = self.global_settings["epsg_code"]
        arr["the_geom"] = pygeos.apply(
            arr["the_geom"], _get_reproject_func(SOURCE_EPSG, target_epsg)
        )

        # transform to a dict of 1D ndarrays
        return ConnectionNodes(**{name: arr[name] for name in arr.dtype.names})

    def get_grid_refinements(self):
        """Return Gridrefinement and GridRefinementArea concatenated into one array.

        """
        with self.get_session() as session:
            data = session.query(
                ST_AsBinary(models.GridRefinement.the_geom),
                models.GridRefinement.display_name,
                models.GridRefinement.id,
                models.GridRefinement.code,
                cast(models.GridRefinement.refinement_level , Integer),
            ).all()

            data += session.query(
                ST_AsBinary(models.GridRefinementArea.the_geom),
                models.GridRefinementArea.display_name,
                models.GridRefinementArea.id,
                models.GridRefinementArea.code,
                cast(models.GridRefinementArea.refinement_level , Integer),

        # transform tuples to a numpy structured array
        arr = np.array(
            data,
            dtype=[
                ("the_geom", "O"),
                ("display_name", "O"),
                ("id", "i8"),
                ("code", "O"),
                ("refinement_level", "i8"),
            ],
        )

        # transform to pygeos.Geometry
        arr["the_geom"] = pygeos.from_wkb(arr["the_geom"])

        # reproject
        target_epsg = self.global_settings["epsg_code"]
        arr["the_geom"] = pygeos.apply(
            arr["the_geom"], _get_reproject_func(SOURCE_EPSG, target_epsg)
        )

        return {name: arr[name] for name in arr.dtype.names}


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
