import geopandas
import pygeos
import numpy as np
import fiona
import pyproj

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
    # use Geopandas for loading (= fiona/GDAL on the background)
    df = geopandas.GeoDataFrame.from_file(path, layer="v2_channel")

    # but go to numpy ASAP to have more control
    geoms = df["geometry"].values.data
    dist_calc_points = df["dist_calc_points"].values

    # reproject
    geoms = pygeos.apply(geoms, reproject)

    # compute number of nodes to add per channel
    length = pygeos.length(geoms) 
    n_segments = np.maximum(
        np.round(length / dist_calc_points).astype(int),
        1,
    )
    segment_size = length(geoms) / n_segments
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


if __name__ == "__main__":
    channel_points(
        "/var/models/v2_bergermeer/3a58dfb2e95f10ca121e26e4259ed68f60e3268f/v2_bergermeer.sqlite",
        out_path="/var/models/noded_channels_pygeos.shp",
    )
