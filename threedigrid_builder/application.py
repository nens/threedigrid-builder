"""The application layer a.k.a. use cases of the project

Use cases orchestrate the flow of data to and from the domain entities.

This layer depends on the interfaces as well as on the domain layer.
"""

from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.grid import Grid
from threedigrid_builder.grid import QuadTree
from threedigrid_builder.interface import SQLite
from threedigrid_builder.interface import Subgrid

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
    grid.set_channel_weights(cross_section_locations, channels)

    grid.epsg_code = db.global_settings["epsg_code"]
    return grid


def get_2d_grid(sqlite_path, dem_path, model_area_path=None):
    """Make 2D computational grid
    """

    subgrid = Subgrid(dem_path, model_area=model_area_path)
    subgrid_meta = subgrid.get_meta()

    db = SQLite(sqlite_path)
    refinements = db.get_grid_refinements()
    quadtree = QuadTree(
        subgrid_meta,
        db.global_settings["kmax"],
        db.global_settings["grid_space"],
        refinements
    )
    grid = Grid.from_quadtree(quadtree)
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

    # construct points from nodes.coordinates
    node_geometries = pygeos.points(node_data.pop("coordinates"))
    # protect against https://github.com/pygeos/pygeos/issues/306
    node_data["bounds"][~np.isfinite(node_data["bounds"])] = 0.
    cell_geometries = pygeos.box(
        node_data["bounds"][:, 0],
        node_data["bounds"][:, 1],
        node_data["bounds"][:, 2],
        node_data["bounds"][:, 3],
    )
    node_data.pop("bounds")

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
    df_cells = geopandas.GeoDataFrame(
        node_data, geometry=cell_geometries, crs=grid.epsg_code
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
    df_cells.to_file(out_path, layer="cells", driver="GPKG")
    df_lines.to_file(out_path, layer="lines", driver="GPKG")
    # TODO make a "line_geometries" layer from lines.line_geometries
