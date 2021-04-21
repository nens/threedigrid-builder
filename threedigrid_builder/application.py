"""The application layer a.k.a. use cases of the project

Use cases orchestrate the flow of data to and from the domain entities.

This layer depends on the interfaces as well as on the domain layer.
"""

from pathlib import Path
from threedigrid_builder.grid import Grid
from threedigrid_builder.grid import QuadTree
from threedigrid_builder.interface import GeopackageOut
from threedigrid_builder.interface import GridAdminOut
from threedigrid_builder.interface import SQLite
from threedigrid_builder.interface import Subgrid
from typing import Optional

import itertools


__all__ = ["make_grid"]


def _get_1d_grid(path, node_id_start=0, line_id_start=0):
    """Compute interpolated channel nodes"""
    db = SQLite(path)

    # the offsets of the node ids are controlled from here
    # for now, we have ConnectionNodes - ChannelNodes:
    connection_nodes = db.get_connection_nodes()

    node_id_counter = itertools.count(start=node_id_start)
    line_id_counter = itertools.count(start=line_id_start)

    grid = Grid.from_connection_nodes(
        connection_nodes=connection_nodes, node_id_counter=node_id_counter
    )

    channels = db.get_channels()
    grid += Grid.from_channels(
        connection_nodes=connection_nodes,
        channels=channels,
        global_dist_calc_points=db.global_settings["dist_calc_points"],
        node_id_counter=node_id_counter,
        line_id_counter=line_id_counter,
        connection_node_offset=grid.nodes.id[0],
    )

    cross_section_locations = db.get_cross_section_locations()
    grid.set_channel_weights(cross_section_locations, channels)
    grid.set_calculation_types()
    grid.set_bottom_levels(cross_section_locations, channels, None, None, None)

    grid.finalize(epsg_code=db.global_settings["epsg_code"])
    return grid


def _get_2d_grid(sqlite_path, dem_path, model_area_path=None):
    """Make 2D computational grid"""

    node_counter = itertools.count()
    line_counter = itertools.count()

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
    grid = Grid.from_quadtree(
        quadtree=quadtree,
        area_mask=subgrid_meta["area_mask"],
        node_id_counter=node_counter,
        line_id_counter=line_counter,
    )

    grid.finalize(epsg_code=db.global_settings["epsg_code"])

    return grid


def _grid_to_gpkg(grid, path):
    with GeopackageOut(path) as out:
        out.write_nodes(grid.nodes, epsg_code=grid.epsg_code)
        out.write_lines(grid.lines, epsg_code=grid.epsg_code)


def _grid_to_hdf5(grid, path):
    with GridAdminOut(path) as out:
        out.write_grid_characteristics(grid.nodes, grid.lines, epsg_code=grid.epsg_code)
        out.write_grid_counts(grid.nodes, grid.lines)
        if grid.quadtree_stats is not None:
            out.write_quadtree(grid.quadtree_stats)
        out.write_nodes(grid.nodes)
        out.write_lines(grid.lines)


def make_grid(
    sqlite_path: Path,
    dem_path: Path,
    out_path: Path,
    model_area_path: Optional[Path] = None,
):
    """Create a Grid instance from sqlite and DEM paths

    The SQLite is expected to be validated already with threedi-modelchecker.

    Args:
        sqlite_path: The path to the input schematisation (SQLite) file
        dem_path: The path of the input DEM file (GeoTIFF)
        out_path: The path of the (to be created) output file. Allowed extensions
            are: .h5 (HDF5) and .gpkg (Geopackage)
        model_area_path

    Raises:
        threedigrid_builder.SchematisationError: if there is something wrong with
            the input schematisation (SQLite) file.
    """
    if isinstance(sqlite_path, str):
        sqlite_path = Path(sqlite_path)
    if isinstance(dem_path, str):
        dem_path = Path(dem_path)
    if isinstance(out_path, str):
        out_path = Path(out_path)
    if isinstance(model_area_path, str):
        model_area_path = Path(model_area_path)
    extension = out_path.suffix.lower()
    if extension == ".h5":
        writer = _grid_to_hdf5
    elif extension:
        writer = _grid_to_gpkg
    else:
        raise ValueError(f"Unsupported output format '{extension}'")

    grid = _get_2d_grid(sqlite_path, dem_path, model_area_path)

    node_id_start = grid.nodes.id[-1] + 1 if len(grid.nodes) > 0 else 0
    line_id_start = grid.lines.id[-1] + 1 if len(grid.lines) > 0 else 0
    grid += _get_1d_grid(
        sqlite_path, node_id_start=node_id_start, line_id_start=line_id_start
    )

    writer(grid, out_path)
