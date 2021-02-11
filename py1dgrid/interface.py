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


@lru_cache(maxsize=8)  # this function is deterministic
def get_reproject_func(source_epsg, target_epsg):
    transformer = Transformer.from_crs(
        CRS.from_epsg(source_epsg), CRS.from_epsg(target_epsg)
    )

    def func(coords):
        x, y = transformer.transform(coords[:, 0], coords[:, 1])
        return np.array([x, y]).T

    return func


def object_as_dict(obj):
    # https://stackoverflow.com/questions/1958219/convert-sqlalchemy-row-object-to-python-dict
    return {c.key: getattr(obj, c.key) for c in inspect(obj).mapper.column_attrs}


@contextmanager
def get_session(path):
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


@lru_cache(maxsize=8)  # the sqlite is assumed not to change!
def get_global_settings(path):
    """Return the global settings from the SQLite at path
    
    Returns:
      dict
    """
    with get_session(path) as session:
        settings = session.query(models.GlobalSetting).one()
    return object_as_dict(settings)


def get_channels(path):
    """Return channels as a dict 1D ndarrays

    Returns:
      dict with the following keys:
      - the_geom (ndarray of pygeos.Geometry)
      - dist_calc_points (ndarray of float32)
      - code (ndarray of python str objects)
      - display_name (ndarray of python str objects)
      - calculation_type (ndarray of uint8)
    """
    with get_session(path) as session:
        data = session.query(
            ST_AsBinary(models.Channel.the_geom),
            models.Channel.dist_calc_points,
            models.Channel.code,
            models.Channel.display_name,
            cast(models.Channel.calculation_type, Integer),
        ).all()

    # transform tuples to a numpy structured array
    arr = np.array(
        data,
        dtype=[
            ("the_geom", "O"),
            ("dist_calc_points", "f4"),
            ("code", "O"),
            ("display_name", "O"),
            ("calculation_type", "u1"),
        ],
    )
    # transform to pygeos.Geometry
    arr["the_geom"] = pygeos.from_wkb(arr["the_geom"])

    # reproject
    target_epsg = get_global_settings(path).epsg_code
    arr["the_geom"] = pygeos.apply(
        arr["the_geom"], get_reproject_func(SOURCE_EPSG, target_epsg)
    )

    # transform to a dict of 1D ndarrays
    return {name: arr[name] for name in arr.dtype.names}
