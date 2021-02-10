from threedi_modelchecker.threedi_database import ThreediDatabase
from threedi_modelchecker.threedi_model import models
from sqlalchemy import cast, Integer
from geoalchemy2.functions import ST_AsBinary
import geopandas
import numpy as np
import pygeos
import fiona


def read_sqlalchemy(path):
    # use SQLAlchemy for loading
    sqlite_settings = {"db_path": path, "db_file": path}
    db = ThreediDatabase(
        connection_settings=sqlite_settings, db_type="spatialite", echo=False
    )
    session = db.get_session()

    data = session.query(
        ST_AsBinary(models.Channel.the_geom),
        models.Channel.dist_calc_points,
        models.Channel.code,
        models.Channel.display_name,
        cast(models.Channel.calculation_type, Integer),
    ).all()

    session.close()

    # transform tuples to a numpy structured array
    arr = np.array(
        data, dtype=[
            ("the_geom", "O"),
            ("dist_calc_points", "f4"),
            ("code", "O"),
            ("display_name", "O"),
            ("calculation_type", "u1"),
        ]
    )

    # transform to dict of 1D arrays for apples-to-apples comparison
    result = {name: arr[name] for name in arr.dtype.names}

    # cast the geometry column
    result["the_geom"] = pygeos.from_wkb(result["the_geom"])

    return result


def read_fiona(path):
    fiona.supported_drivers["SQLite"] = "r"
    # geopandas has an optimized read implementation
    df = geopandas.GeoDataFrame.from_file(
        path,
        layer="v2_channel",
        ignore_fields=[
            "zoom_category",
            "connection_node_start_id",
            "connection_node_end_id",
        ]
    )

    return {
        "the_geom": df.geometry.values.data,
        "dist_calc_points": df["dist_calc_points"].values,
        "code": df["code"].values,
        "display_name": df["display_name"].values,
        "calculation_type": df["calculation_type"].values,
    }
