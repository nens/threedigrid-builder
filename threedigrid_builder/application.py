from threedigrid_builder.interface import SQLite


def export_1D_nodes_to_shapefile(path, out_path):
    """Compute interpolated channel nodes and export to shapefile"""
    import geopandas

    db = SQLite(path)

    channels = db.get_channels()
    points = channels.interpolate_nodes(db.global_settings["dist_calc_points"])

    df = geopandas.GeoDataFrame(points, crs=db.global_settings["epsg_code"])
    df.set_geometry("geometry")
    df.to_file(out_path)
