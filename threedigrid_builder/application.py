from threedigrid_builder.interface import SQLite


def get_1d_grid(path):
    """Compute interpolated channel nodes"""
    db = SQLite(path)

    connection_nodes = db.get_connection_nodes()
    channels = db.get_channels()

    channel_nodes = channels.interpolate_nodes(
        global_dist_calc_points=db.global_settings["dist_calc_points"]
    )

    connection_node_grid = connection_nodes.get_grid()
    channel_grid = channels.get_grid(
        channel_nodes,
        connection_node_offset=0,
        channel_node_offset=len(connection_node_grid.nodes),
    )

    grid = connection_node_grid.concatenate(channel_grid)

    # import geopandas
    # df_nodes = geopandas.GeoDataFrame(
    #     geometry=pygeos.points()

    #     , crs=db.global_settings["epsg_code"])
    # df.set_geometry("geometry")
    # df.to_file(out_path)
