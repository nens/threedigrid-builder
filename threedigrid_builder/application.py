from threedigrid_builder.grid1d import node_channels
from threedigrid_builder.interface import get_global_settings


def export_1D_nodes_to_shapefile(path, out_path):
    """Compute interpolated channel nodes and export to shapefile"""

    import geopandas

    points = node_channels(path=path)
    df = geopandas.GeoDataFrame(points, crs=get_global_settings(path)["epsg_code"])
    df.set_geometry("geometry")
    df.to_file(out_path)
