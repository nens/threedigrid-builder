"""The application layer a.k.a. use cases of the project

Use cases orchestrate the flow of data to and from the domain entities.

This layer depends on the interfaces as well as on the domain layer.
"""

from pathlib import Path
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import Grid
from threedigrid_builder.grid import QuadTree
from threedigrid_builder.interface import GeopackageOut
from threedigrid_builder.interface import GridAdminOut
from threedigrid_builder.interface import SQLite
from threedigrid_builder.interface import Subgrid
from typing import Callable
from typing import Optional

import itertools
import logging


__all__ = ["make_grid"]

logger = logging.getLogger(__name__)


def _default_progress_callback(progress: float, message: str):
    logger.info("Progress: %d, Message: %s", progress * 100, message)


def _make_grid(
    sqlite_path,
    dem_path=None,
    model_area_path=None,
    meta=None,
    progress_callback=None,
):
    """Compute interpolated channel nodes"""
    progress_callback(0.0, "Reading input schematisation...")
    db = SQLite(sqlite_path)

    node_id_counter = itertools.count()
    line_id_counter = itertools.count()

    settings = db.get_settings()
    grid = Grid.from_meta(**settings, **(meta or {}))
    grid_settings = settings["grid_settings"]

    if grid_settings.use_2d:
        progress_callback(0.1, "Constructing subgrid...")
        if not dem_path:
            raise SchematisationError("DEM file expected")
        # TODO use_2d_flow --> https://github.com/nens/threedigrid-builder/issues/87
        subgrid = Subgrid(dem_path, model_area=model_area_path)
        subgrid_meta = subgrid.get_meta()
        refinements = db.get_grid_refinements()
        progress_callback(0.7, "Constructing quadtree...")
        quadtree = QuadTree(
            subgrid_meta,
            grid_settings.kmax,
            grid_settings.grid_space,
            refinements,
        )
        grid += Grid.from_quadtree(
            quadtree=quadtree,
            area_mask=subgrid_meta["area_mask"],
            node_id_counter=node_id_counter,
            line_id_counter=line_id_counter,
        )

    connection_nodes = db.get_connection_nodes()
    if grid_settings.use_1d_flow and len(connection_nodes) > 0:
        progress_callback(0.8, "Constructing 1D network...")
        cn_grid = Grid.from_connection_nodes(
            connection_nodes=connection_nodes, node_id_counter=node_id_counter
        )
        connection_node_first_id = cn_grid.nodes.id[0] if len(cn_grid.nodes) > 0 else 0
        grid += cn_grid

        channels = db.get_channels()
        grid += Grid.from_channels(
            connection_nodes=connection_nodes,
            channels=channels,
            global_dist_calc_points=grid_settings.dist_calc_points,
            node_id_counter=node_id_counter,
            line_id_counter=line_id_counter,
            connection_node_offset=connection_node_first_id,
        )

        cross_section_locations = db.get_cross_section_locations()
        grid.set_channel_weights(cross_section_locations, channels)

        pipes = db.get_pipes()
        grid += Grid.from_pipes(
            connection_nodes=connection_nodes,
            pipes=pipes,
            global_dist_calc_points=grid_settings.dist_calc_points,
            node_id_counter=node_id_counter,
            line_id_counter=line_id_counter,
            connection_node_offset=connection_node_first_id,
        )

        culverts = db.get_culverts()
        weirs = db.get_weirs()
        orifices = db.get_orifices()
        grid += grid.from_structures(
            connection_nodes=connection_nodes,
            culverts=culverts,
            weirs=weirs,
            orifices=orifices,
            global_dist_calc_points=grid_settings.dist_calc_points,
            node_id_counter=node_id_counter,
            line_id_counter=line_id_counter,
            connection_node_offset=connection_node_first_id,
        )

        grid.set_calculation_types()
        grid.set_bottom_levels(
            cross_section_locations, channels, pipes, weirs, orifices, culverts
        )
        grid.set_pumps(db.get_pumps())

    if grid.meta.has_1d and grid.meta.has_2d:
        progress_callback(0.9, "Connecting 1D and 2D elements...")
        grid.add_1d2d(
            connection_nodes=connection_nodes,
            channels=channels,
            pipes=pipes,
            locations=cross_section_locations,
            culverts=culverts,
            line_id_counter=line_id_counter,
        )

    grid.finalize()
    return grid


def _grid_to_gpkg(grid, path):
    with GeopackageOut(path) as out:
        out.write_nodes(grid.nodes, epsg_code=grid.meta.epsg_code)
        out.write_lines(grid.lines, epsg_code=grid.meta.epsg_code)
        out.write_pumps(grid.pumps, epsg_code=grid.meta.epsg_code)


def _grid_to_hdf5(grid, path):
    with GridAdminOut(path) as out:
        out.write_meta(grid.meta)
        out.write_grid_counts(grid.nodes, grid.lines)
        if grid.quadtree_stats is not None:
            out.write_quadtree(grid.quadtree_stats)
        out.write_nodes(grid.nodes)
        out.write_lines(grid.lines)
        out.write_pumps(grid.pumps)


def make_grid(
    sqlite_path: Path,
    dem_path: Path,
    out_path: Path,
    model_area_path: Optional[Path] = None,
    meta: dict = None,
    progress_callback: Optional[Callable[[float, str], None]] = None,
):
    """Create a Grid instance from sqlite and DEM paths

    The SQLite is expected to be validated already with threedi-modelchecker.

    Args:
        sqlite_path: The path to the input schematisation (SQLite) file
        dem_path: The path of the input DEM file (GeoTIFF)
        out_path: The path of the (to be created) output file. Allowed extensions
            are: .h5 (HDF5) and .gpkg (Geopackage)
        model_area_path
        meta: an optional dict with the following (optional) keys: model_slug (str),
            revision_hash (str), revision_nr (int), threedi_version (str)
        progress_callback: an optional function that updates the progress. The function
            should take an float in the range 0-1 and a message string.

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
    if progress_callback is None:
        progress_callback = _default_progress_callback

    grid = _make_grid(
        sqlite_path,
        dem_path,
        model_area_path,
        meta=meta,
        progress_callback=progress_callback,
    )

    progress_callback(0.95, "Writing gridadmin...")
    writer(grid, out_path)
