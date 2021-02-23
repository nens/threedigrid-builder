# Module containing functions to read from the SQLite

from condenser import NumpyQuery
from contextlib import contextmanager
from functools import lru_cache
from geoalchemy2.functions import ST_AsBinary
from pyproj import Transformer
from pyproj.crs import CRS
from sqlalchemy.orm import sessionmaker
from sqlalchemy import cast
from sqlalchemy import inspect
from sqlalchemy import Float
from sqlalchemy import Integer
from threedi_modelchecker.threedi_database import ThreediDatabase
from threedi_modelchecker.threedi_model import models
from threedigrid_builder.grid1d import Channels
from threedigrid_builder.grid1d import ConnectionNodes

import numpy as np
import pygeos


from threedi_modelchecker.threedi_model.custom_types import IntegerEnum


SOURCE_EPSG = 4326

__all__ = ["SQLite"]

# put some global defaults on datatypes
NumpyQuery.default_numpy_settings[Float]["dtype"] = np.float32
NumpyQuery.default_numpy_settings[Integer]["dtype"] = np.int32
NumpyQuery.default_numpy_settings[IntegerEnum] = {
    "sql_cast": lambda x: cast(x, Integer),
    "dtype": np.int32,
}


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
        session = self.db.get_session(query_cls=NumpyQuery)
        try:
            yield session
        finally:
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
        """
        target_epsg = self.global_settings["epsg_code"]
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
                .with_transformed_geometries(target_epsg)
                .as_structarray()
            )

        # reproject
        # target_epsg = self.global_settings["epsg_code"]
        # arr["the_geom"] = pygeos.apply(
        #     arr["the_geom"], _get_reproject_func(SOURCE_EPSG, target_epsg)
        # )

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
