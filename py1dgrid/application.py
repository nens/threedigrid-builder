from .domain import get_target_epsg, node_channels


def export_1D_nodes_to_shapefile(path, out_path):
    """Compute additional channel nodes and export to shapefile"""

    import geopandas

    points = node_channels(path=path)
    df = geopandas.GeoDataFrame(points, crs=get_target_epsg(path))
    df.set_geometry("geometry")
    df.to_file(out_path)
