from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.grid import Grid
from threedigrid_builder.interface import SQLite

import itertools
import numpy as np
import pygeos


def get_1d_grid(path):
    """Compute interpolated channel nodes"""
    db = SQLite(path)

    # the offsets of the node ids are controlled from here
    # for now, we have ConnectionNodes - ChannelNodes:
    connection_nodes = db.get_connection_nodes()
    counter = itertools.count()

    grid = Grid.from_connection_nodes(
        connection_nodes=connection_nodes, node_id_counter=counter
    )

    channels = db.get_channels()
    grid += Grid.from_channels(
        connection_nodes=connection_nodes,
        channels=channels,
        global_dist_calc_points=db.global_settings["dist_calc_points"],
        node_id_counter=counter,
    )

    cross_section_locations = db.get_cross_section_locations()
    cross_section_locations.apply_to_channels(channels, grid.lines)

    grid.epsg_code = db.global_settings["epsg_code"]
    return grid


def _enum_to_str(arr, enum_type):
    result = np.full_like(arr, "", dtype=object)
    result[arr != -9999] = [enum_type(x).name for x in arr[arr != -9999]]
    return result


def grid_to_gpkg(grid, out_path):
    """Export grid to a geopackage. Geopandas is required."""
    import geopandas  # optional dependency for easy file output

    # ---NODES---

    node_data = grid.nodes.to_dict()
    node_data.pop("bounds")  # cannot export 2D arrays

    # construct points from nodes.coordinates
    node_geometries = pygeos.points(node_data.pop("coordinates"))

    # convert enums to strings
    node_data["node_type"] = _enum_to_str(node_data["node_type"], NodeType)
    node_data["content_type"] = _enum_to_str(node_data["content_type"], ContentType)
    node_data["calculation_type"] = _enum_to_str(
        node_data["calculation_type"], CalculationType
    )

    # construct the geodataframe
    df_nodes = geopandas.GeoDataFrame(
        node_data, geometry=node_geometries, crs=grid.epsg_code
    )

    # ---LINES---

    line_data = grid.lines.to_dict()
    line_data.pop("line_geometries")  # cannot export geometries

    # construct linestrings from nodes.coordinates and lines.line
    line_ids = line_data.pop("line")
    start = grid.nodes.coordinates[line_ids[:, 0]]
    end = grid.nodes.coordinates[line_ids[:, 1]]
    line_geometries = pygeos.linestrings(
        np.concatenate([start[:, np.newaxis], end[:, np.newaxis]], axis=1)
    )

    # convert enums to strings
    line_data["line_type"] = _enum_to_str(line_data["line_type"], LineType)
    line_data["content_type"] = _enum_to_str(line_data["content_type"], ContentType)

    # construct the geodataframe
    df_lines = geopandas.GeoDataFrame(
        line_data, geometry=line_geometries, crs=grid.epsg_code
    )

    # ---WRITE THE FILE---
    df_nodes.to_file(out_path, layer="nodes", driver="GPKG")
    df_lines.to_file(out_path, layer="lines", driver="GPKG")
    # TODO make a "cells" layer from nodes.bounds
    # TODO make a "line_geometries" layer from lines.line_geometries
