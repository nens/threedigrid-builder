"""The application layer a.k.a. use cases of the project

Use cases orchestrate the flow of data to and from the domain entities.

This layer depends on the interfaces as well as on the domain layer.
"""

from threedigrid_builder.grid import Grid
from threedigrid_builder.grid import QuadTree
from threedigrid_builder.interface import GeopackageOut
from threedigrid_builder.interface import GridAdminOut
from threedigrid_builder.interface import SQLite
from threedigrid_builder.interface import Subgrid

import itertools


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

    grid.finalize(epsg_code=db.global_settings["epsg_code"], pixel_size=None)
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
        refinements,
    )
    grid = Grid.from_quadtree(quadtree, subgrid_meta)

    grid.finalize(
        epsg_code=db.global_settings["epsg_code"], pixel_size=subgrid_meta["pixel_size"]
    )
    return grid


def grid_to_gpkg(grid, path):
    with GeopackageOut(path) as out:
        out.write_nodes(grid.nodes, epsg_code=grid.epsg_code)
        out.write_lines(grid.lines, epsg_code=grid.epsg_code)


def grid_to_hdf5(grid, path):
    with GridAdminOut(path) as out:
        out.write_nodes(grid.nodes, pixel_size=grid.pixel_size)
        out.write_lines(grid.lines, epsg_code=grid.epsg_code)
