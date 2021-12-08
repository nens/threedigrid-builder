"""The application layer a.k.a. use cases of the project

Use cases orchestrate the flow of data to and from the domain entities.

This layer depends on the interfaces as well as on the domain layer.
"""

from pathlib import Path
from threedigrid_builder.base.surfaces import Surfaces
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import Grid
from threedigrid_builder.grid import QuadTree
from threedigrid_builder.grid.zero_d import ImperviousSurfaces
from threedigrid_builder.interface import DictOut
from threedigrid_builder.interface import GDALInterface
from threedigrid_builder.interface import GeopackageOut
from threedigrid_builder.interface import GridAdminOut
from threedigrid_builder.interface import RasterioInterface
from threedigrid_builder.interface import SQLite
from threedigrid_builder.constants import InflowType
from typing import Callable
from typing import Optional

import itertools
import logging


__all__ = ["make_grid", "make_gridadmin"]

logger = logging.getLogger(__name__)


def _default_progress_callback(progress: float, message: str):
    logger.info("Progress: %d, Message: %s", progress * 100, message)


def _make_gridadmin(
    sqlite_path,
    dem_path=None,
    meta=None,
    progress_callback=None,
):
    """Compute interpolated channel nodes"""
    progress_callback(0.0, "Reading input schematisation...")
    db = SQLite(sqlite_path)

    node_id_counter = itertools.count()
    line_id_counter = itertools.count()
    embedded_node_id_counter = itertools.count()

    settings = db.get_settings()
    grid = Grid.from_meta(**settings, **(meta or {}))
    grid_settings = settings["grid_settings"]

    if grid_settings.use_2d:
        progress_callback(0.1, "Reading subgrid input...")
        if not dem_path:
            raise SchematisationError("DEM file expected")
        # TODO use_2d_flow --> https://github.com/nens/threedigrid-builder/issues/87

        try:
            with RasterioInterface(dem_path) as raster:
                subgrid_meta = raster.read()
        except ImportError:
            with GDALInterface(dem_path) as raster:
                subgrid_meta = raster.read()

        # Patch epsg code with that of the DEM (so: user-supplied EPSG is ignored)
        grid.meta.epsg_code = raster.epsg_code

        refinements = db.get_grid_refinements()
        progress_callback(0.7, "Constructing 2D computational grid...")
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
        grid.set_obstacles(db.get_obstacles(), db.get_levees())
        grid.set_boundary_conditions_2d(
            db.get_boundary_conditions_2d(),
            quadtree,
            node_id_counter,
            line_id_counter,
        )
        dem_average_areas = db.get_dem_average_areas()
        grid.set_dem_averaged_cells(dem_average_areas)

    connection_nodes = db.get_connection_nodes()
    if grid_settings.use_1d_flow and len(connection_nodes) > 0:
        progress_callback(0.8, "Constructing 1D computational grid...")
        cn_grid = Grid.from_connection_nodes(
            connection_nodes=connection_nodes,
            node_id_counter=node_id_counter,
        )
        connection_node_first_id = cn_grid.nodes.id[0] if len(cn_grid.nodes) > 0 else 0
        grid += cn_grid

        channels = db.get_channels()
        channel_grid = Grid.from_linear_objects(
            connection_nodes=connection_nodes,
            objects=channels,
            cell_tree=grid.cell_tree if grid_settings.use_2d else None,
            global_dist_calc_points=grid_settings.dist_calc_points,
            embedded_cutoff_threshold=grid_settings.embedded_cutoff_threshold,
            node_id_counter=node_id_counter,
            embedded_node_id_counter=embedded_node_id_counter,
            line_id_counter=line_id_counter,
            connection_node_offset=connection_node_first_id,
        )

        locations = db.get_cross_section_locations()
        locations.apply_to_lines(channel_grid.lines, channels)
        grid += channel_grid

        pipes = db.get_pipes()
        grid += Grid.from_linear_objects(
            connection_nodes=connection_nodes,
            objects=pipes,
            cell_tree=grid.cell_tree if grid_settings.use_2d else None,
            global_dist_calc_points=grid_settings.dist_calc_points,
            embedded_cutoff_threshold=grid_settings.embedded_cutoff_threshold,
            node_id_counter=node_id_counter,
            embedded_node_id_counter=embedded_node_id_counter,
            line_id_counter=line_id_counter,
            connection_node_offset=connection_node_first_id,
        )

        culverts = db.get_culverts()
        grid += Grid.from_linear_objects(
            connection_nodes=connection_nodes,
            objects=culverts,
            cell_tree=grid.cell_tree if grid_settings.use_2d else None,
            global_dist_calc_points=grid_settings.dist_calc_points,
            embedded_cutoff_threshold=grid_settings.embedded_cutoff_threshold,
            node_id_counter=node_id_counter,
            embedded_node_id_counter=embedded_node_id_counter,
            line_id_counter=line_id_counter,
            connection_node_offset=connection_node_first_id,
        )

        weirs = db.get_weirs()
        orifices = db.get_orifices()
        grid += grid.from_structures(
            connection_nodes=connection_nodes,
            weirs=weirs,
            orifices=orifices,
            line_id_counter=line_id_counter,
            connection_node_offset=connection_node_first_id,
        )

        grid.set_calculation_types()
        grid.set_bottom_levels()
        grid.set_initial_waterlevels(
            connection_nodes=connection_nodes,
            channels=channels,
            pipes=pipes,
            culverts=culverts,
        )
        grid.set_boundary_conditions_1d(db.get_boundary_conditions_1d())
        grid.set_cross_sections(db.get_cross_section_definitions())
        grid.set_pumps(db.get_pumps())

    if grid.nodes.has_1d and grid.nodes.has_2d:
        progress_callback(0.9, "Connecting 1D and 2D domains...")
        grid.embed_nodes(embedded_node_id_counter)
        connected_points = db.get_connected_points()
        grid.add_1d2d(
            connected_points,
            connection_nodes=connection_nodes,
            channels=channels,
            pipes=pipes,
            locations=locations,
            culverts=culverts,
            line_id_counter=line_id_counter,
        )
        grid.add_breaches(connected_points)

    if grid_settings.use_0d_inflow in (
        InflowType.IMPERVIOUS_SURFACE.value,
        InflowType.SURFACE.value,
    ):
        # process zero-d information, either using the 'v2_surfaces'
        # or 'v2_impervious_surfaces' table.
        progress_callback(0.95, "Processing 0D domain...")
        if grid_settings.use_0d_inflow == InflowType.SURFACE.value:
            surfaces: Surfaces = db.get_surfaces()
        else:
            surfaces: ImperviousSurfaces = db.get_impervious_surfaces()
        grid.add_0d(surfaces)

    grid.finalize()
    return grid


def make_gridadmin(
    sqlite_path: Path,
    dem_path: Path,
    out_path: Optional[Path] = None,
    meta: dict = None,
    progress_callback: Optional[Callable[[float, str], None]] = None,
):
    """Create a Grid instance from sqlite and DEM paths

    The SQLite is expected to be validated already with threedi-modelchecker.

    Args:
        sqlite_path: The path to the input schematisation (SQLite) file
        dem_path: The path of the input DEM file (GeoTIFF)
        out_path: The path of the (to be created) output file. Allowed extensions
            are: .h5 (HDF5) and .gpkg (Geopackage). If not supplied, this function will
            return an in-memory representation of the Grid.
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
    if out_path is None:
        Writer = DictOut
    else:
        if isinstance(out_path, str):
            out_path = Path(out_path)
        extension = out_path.suffix.lower()
        if extension == ".h5":
            Writer = GridAdminOut
        elif extension == ".gpkg":
            Writer = GeopackageOut
        else:
            raise ValueError(f"Unsupported output format '{extension}'")

    if not Writer.available():
        raise ImportError(
            "This output format is not available because dependencies are missing."
        )
    if progress_callback is None:
        progress_callback = _default_progress_callback

    grid = _make_gridadmin(
        sqlite_path,
        dem_path,
        meta=meta,
        progress_callback=progress_callback,
    )

    progress_callback(0.99, "Writing gridadmin...")
    with Writer(out_path) as writer:
        return writer.write(grid)


# Legacy
make_grid = make_gridadmin
