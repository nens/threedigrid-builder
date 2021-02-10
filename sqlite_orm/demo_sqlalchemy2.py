import geopandas
import pygeos
import numpy as np
import fiona
import pyproj
from sqlalchemy import cast, Integer
from geoalchemy2.functions import ST_AsBinary
from threedi_modelchecker.threedi_database import ThreediDatabase
from threedi_modelchecker.threedi_model import models

# https://github.com/Toblerity/Fiona/issues/992
fiona.supported_drivers["SQLite"] = "r"

# Note: pygeos and shapely will be merged in the course of 2021, until then 
# install from source to match GEOS versions and so fast interoperability
# between geopandas/shapely and pygeos.
#
# pip install pygeos --no-binary pygeos
# pip install shapely --no-binary shapely

transformer = pyproj.Transformer.from_crs(
    "epsg:4326", "epsg:28992", always_xy=True
)

def reproject(coords):
    x, y = transformer.transform(coords[:, 0], coords[:, 1])
    return np.array([x, y]).T


def channel_points(path, out_path=None):
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
    # transform to geometries
    geoms = pygeos.from_wkb(arr["the_geom"])

    # reproject
    geoms = pygeos.apply(geoms, reproject)

    # compute number of nodes to add per channel
    length = pygeos.length(geoms) 
    n_segments = np.maximum(
        np.round(length / arr["dist_calc_points"]).astype(int),
        1,
    )
    segment_size = length / n_segments
    n_nodes = n_segments - 1
    idx = np.repeat(np.arange(geoms.size), n_nodes)

    # some numpy juggling to get the distance to the start of each channel
    dist_to_start = np.arange(idx.size)
    dist_to_start[n_nodes[0]:] -= np.repeat(np.cumsum(n_nodes)[:-1], n_nodes[1:])
    dist_to_start = (dist_to_start + 1) * segment_size[idx]

    points = pygeos.line_interpolate_point(
        geoms[idx],  # note: this only copies geometry pointers
        dist_to_start,
    )

    out = geopandas.GeoDataFrame({
        "geometry": points,
        "ch_code": arr["code"][idx],
        "ch_name": arr["display_name"][idx],
        "ch_type": arr["calculation_type"][idx],
    })

    # transform back to geopandas for IO
    if out_path:
        out.to_file(out_path)

    return out


def channel_conn_nodes(path, projection="epsg:28992", out_path=None):
    df = geopandas.GeoDataFrame.from_file(path, layer="v2_connection_nodes")
    # v2_connection_nodes has 2 different geometry columns, with 2 different
    # 'layers' in the geometry_columns table
    # SELECT * FROM geometry_columns WHERE f_table_name='v2_connection_nodes';
    # [('v2_connection_nodes', 'the_geom', 'POINT', 'XY', 4326, 1),
    #  ('v2_connection_nodes', 'the_geom_linestring', 'LINESTRING', 'XY', 4326, 0)]
    #
    # There seems to be no way to access the points through Fiona! Tried:
    # - accessing by layer name -> yields geometry=None everywhere and in the
    #   schema "Geometry": "LineString"
    # - accessing by layer id -> there is only one with name=v2_connection_nodes
    # - adding various options (propagated to GDALOpenEx) - no option seems
    #   to have an effect


if __name__ == "__main__":
    channel_points(
        "/var/models/v2_bergermeer/3a58dfb2e95f10ca121e26e4259ed68f60e3268f/v2_bergermeer.sqlite",
        out_path="/var/models/noded_channels_pygeos2.shp",
    )
