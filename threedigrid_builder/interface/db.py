# Module containing functions to read from the SQLite

from contextlib import contextmanager
from functools import lru_cache

import numpy as np
import pygeos
from geoalchemy2.functions import ST_AsBinary
from pyproj import Transformer
from pyproj.crs import CRS
from sqlalchemy import Integer, cast, inspect
from threedi_modelchecker.threedi_database import ThreediDatabase
from threedi_modelchecker.threedi_model import models

SOURCE_EPSG = 4326

__all__ = ["get_global_settings", "get_channels"]


@contextmanager
def _get_session(path):
    """A context manager that yields an SQLAlchemy session.

    The session is closed af the context manager exit. No commit or rollback
    is done. It is meant for read-only access, but writing is not prohibited.

    Args:
      path (str): Path to an SQLite
    
    Returns:
      SQLAlchemy.orm.Session
    """
    sqlite_settings = {"db_path": path, "db_file": path}
    db = ThreediDatabase(
        connection_settings=sqlite_settings, db_type="spatialite", echo=False
    )
    session = db.get_session()
    yield session
    session.close()


def _object_as_dict(obj):
    # https://stackoverflow.com/questions/1958219/convert-sqlalchemy-row-object-to-python-dict
    return {c.key: getattr(obj, c.key) for c in inspect(obj).mapper.column_attrs}


def get_global_settings(path):
    """Return the global settings dictionary from the SQLite at path
    
    Args:
      path (str): Path to an SQLite

    Returns:
      dict
    """
    with _get_session(path) as session:
        settings = session.query(models.GlobalSetting).one()
    return _object_as_dict(settings)


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


def get_channels(path):
    """Return channels as a dict of 1D ndarrays

    Args:
      path (str): Path to an SQLite

    Returns:
      dict with the following keys:
      - the_geom (ndarray of pygeos.Geometry)
      - dist_calc_points (ndarray of float32)
      - code (ndarray of python str objects)
      - display_name (ndarray of python str objects)
      - calculation_type (ndarray of uint8)

    Note:
      Most stuff in this function will eventually be implemented in the package
      "condenser".
    """
    with _get_session(path) as session:
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
    target_epsg = get_global_settings(path)["epsg_code"]
    arr["the_geom"] = pygeos.apply(
        arr["the_geom"], _get_reproject_func(SOURCE_EPSG, target_epsg)
    )

    # transform to a dict of 1D ndarrays
    return {name: arr[name] for name in arr.dtype.names}
