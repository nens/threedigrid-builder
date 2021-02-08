import geopandas
import pygeos
import numpy as np
import fiona

# https://github.com/Toblerity/Fiona/issues/992
fiona.supported_drivers["SQLite"] = "r"

# Note: pygeos and shapely will be merged in the course of 2021, until then 
# install from source to match GEOS versions and so fast interoperability
# between geopandas/shapely and pygeos.
#
# pip install pygeos --no-binary pygeos
# pip install shapely --no-binary shapely

def channel_points(path, projection="epsg:28992", out_path=None):
    # use Geopandas for loading (= fiona/GDAL on the background)
    df = geopandas.GeoDataFrame.from_file(path, layer="v2_channel")

    # but go to numpy ASAP to have more control
    geoms = df["geometry"].to_crs(projection).values.data
    dist_calc_points = df["dist_calc_points"].values

    # compute number of nodes to add per channel
    # TODO This does not give correct results. Previous version probably has
    # some more logic than this (e.g. I get exactly 100m between nodes while
    # previous version seems to distribute them evenly on the channel)
    n_nodes = np.floor(pygeos.length(geoms) / dist_calc_points).astype(int)
    idx = np.repeat(np.arange(geoms.size), n_nodes)

    # per node, compute the distance to the start of the channel 
    # TODO vectorize this using
    # https://stackoverflow.com/questions/49178977/multiple-cumulative-sum-within-a-numpy-array
    dist_to_start = np.empty(n_nodes.sum(), dtype=float)
    cur = 0
    for count, dist_between in zip(n_nodes, dist_calc_points):
        dist_to_start[cur:cur+count] = [(i + 1) * dist_between for i in range(count)]
        cur += count

    points = pygeos.line_interpolate_point(
        geoms[idx],  # note: this only copies geometry pointers
        dist_to_start,
    )

    out = geopandas.GeoDataFrame({
        "geometry": points,
        "ch_code": df.loc[idx, "code"],
        "ch_name": df.loc[idx, "display_name"],
        "ch_type": df.loc[idx, "calculation_type"],
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
