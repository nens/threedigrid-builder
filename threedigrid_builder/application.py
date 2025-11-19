"""The application layer a.k.a. use cases of the project

Use cases orchestrate the flow of data to and from the domain entities.

This layer depends on the interfaces as well as on the domain layer.
"""
import itertools
import logging
from pathlib import Path
from typing import Callable, Optional

import numpy as np
from osgeo import gdal

from threedigrid_builder.base import GridSettings
from threedigrid_builder.base.surfaces import Surfaces
from threedigrid_builder.constants import InflowType, NO_DATA_VALUE
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import Clone, Grid, QuadTree
from threedigrid_builder.grid.fragments import Fragments
from threedigrid_builder.interface import (
    DictOut,
    GDALInterface,
    GeopackageOut,
    GridAdminOut,
    SQLite,
)
from threedigrid_builder.interface.db import map_cross_section_definition

__all__ = ["make_grid", "make_gridadmin"]

logger = logging.getLogger(__name__)


def _default_progress_callback(progress: float, message: str):
    logger.info("Progress: %d, Message: %s", progress * 100, message)


def _make_gridadmin(
    sqlite_path,
    dem_path=None,
    meta=None,
    progress_callback=None,
    upgrade=False,
    convert_to_geopackage=False,
    apply_cutlines=False,
):
    """Compute interpolated channel nodes"""
    progress_callback(0.0, "Reading input schematisation...")
    db = SQLite(
        sqlite_path, upgrade=upgrade, convert_to_geopackage=convert_to_geopackage
    )

    node_id_counter = itertools.count()
    line_id_counter = itertools.count()
    embedded_node_id_counter = itertools.count()

    settings = db.get_settings()
    grid = Grid.from_meta(**settings, **(meta or {}))
    grid_settings: GridSettings = settings["grid_settings"]

    if grid_settings.use_2d:
        progress_callback(0.1, "Reading subgrid input...")
        if not dem_path:
            raise SchematisationError("DEM file expected")
        with GDALInterface(dem_path) as raster:
            subgrid_meta = raster.read()

        # Patch CRS with that of the DEM (so: user-supplied EPSG is ignored)
        grid.set_crs(raster.crs)
        db.epsg_code = grid.meta.epsg_code

        refinements = db.get_grid_refinements()
        progress_callback(0.7, "Constructing 2D computational grid...")
        quadtree = QuadTree(
            subgrid_meta,
            grid_settings.kmax,
            grid_settings.grid_space,  # min gridsize in meters
            grid_settings.use_2d_flow,
            refinements,
        )

        grid += Grid.from_quadtree(
            quadtree=quadtree,
            area_mask=subgrid_meta["area_mask"],
            node_id_counter=node_id_counter,
            line_id_counter=line_id_counter,
        )

        if apply_cutlines:
            (
                fragment_mask,
                node_fragment_array,
                fragment_geometries,
            ) = grid.apply_cutlines(db.get_obstacles(), dem_path)
            # Flip and transpose mask to mimic GDALInterface.read()
            fortran_fragment_mask = np.asfortranarray(np.flipud(fragment_mask).T)
            fortran_node_fragment_array = np.asfortranarray(node_fragment_array)

            size = len(fragment_geometries)
            fragments_centroid = np.empty((size, 2), dtype=np.float64, order="F")
            fragments_polygon0 = []
            max_corners = 0
            for frg_idx in fragment_geometries:
                fragments_centroid[frg_idx] = np.array(
                    fragment_geometries[frg_idx].centroid.coords
                )
                fragments_polygon0.append(fragment_geometries[frg_idx].boundary.coords)
                max_corners = max(max_corners, len(fragments_polygon0[frg_idx].xy[0]))

            fragments_polygon1 = []
            for frg_idx in fragment_geometries:
                pad_width = max_corners - int(len(fragments_polygon0[frg_idx].xy[0]))
                a = np.array(fragments_polygon0[frg_idx])
                padded_array = np.pad(
                    a, ((0, pad_width), (0, 0)), "constant", constant_values=0.0
                )
                fragments_polygon1.append(padded_array)

            fragments_polygon = np.array(fragments_polygon1, dtype=np.float64)

            # Regenerate the nodes and lines including clone cells
            clone = Clone(
                fortran_node_fragment_array,
                fortran_fragment_mask,
                fragments_centroid,
                fragments_polygon,
                quadtree,
                grid,
                area_mask=subgrid_meta["area_mask"],
            )

            # Reset the node/line counter
            node_id_counter = itertools.count()
            line_id_counter = itertools.count()

            # Overwrite nodes and lines attributes according to the new cells and lines numbers
            grid.nodes, grid.lines = Grid.from_clone(
                quadtree,
                clone,
                node_id_counter,
                line_id_counter,
            )
            # Renumber fragment mask and node_fragment_array and store Fragments to model for export
            # The gridadmin/gpkg are exported 1-based (Fortran), also export mask this way (update mapping)
            clone_mapping = clone.clone_numbering.copy()
            for idx in range(len(clone_mapping)):
                if clone_mapping[idx] != NO_DATA_VALUE:
                    clone_mapping[idx] += 1

            original_fragment_mask = (
                fragment_mask.copy()
            )  # We need a copy to prevent overwrite
            for old_fragment_idx, new_fragment_idx in enumerate(clone_mapping):
                fragment_mask[
                    original_fragment_mask == old_fragment_idx
                ] = new_fragment_idx
            # Export fragment tiff
            with GDALInterface(dem_path) as raster:
                subgrid_meta = raster.read()
                fragment_path = (
                    # dem_path.parent / f"fragments_{dem_path.parent.parent.name}.tif"
                    dem_path.parent
                    / "fragments.tif"
                )
                target_ds = gdal.GetDriverByName("GTiff").Create(
                    str(fragment_path),
                    subgrid_meta["width"],
                    subgrid_meta["height"],
                    1,
                    gdal.GDT_Int32,
                )
                target_ds.SetGeoTransform(raster._dataset.GetGeoTransform())
                target_ds.SetProjection(raster._dataset.GetProjection())
                target_ds.GetRasterBand(1).SetNoDataValue(NO_DATA_VALUE)
                target_ds.GetRasterBand(1).WriteArray(fragment_mask)

            # Export to model (and h5/gpkg)
            fragment_ids = []
            fragment_geoms = []
            for n in range(node_fragment_array.shape[0]):
                for f in range(node_fragment_array.shape[1]):
                    fragment_id = node_fragment_array[n][f]
                    if fragment_id != NO_DATA_VALUE:
                        fragment_geom = fragment_geometries[fragment_id]
                        # Map the final index
                        mapped_fragment_id = clone.clone_numbering[fragment_id]
                        if mapped_fragment_id != NO_DATA_VALUE:
                            fragment_ids.append(mapped_fragment_id)
                            fragment_geoms.append(fragment_geom)
            grid.fragments = Fragments(id=fragment_ids, the_geom=fragment_geoms)

        grid.set_quarter_administration(quadtree)

        if grid.meta.has_groundwater:
            grid.add_groundwater(
                grid.meta.has_groundwater_flow, node_id_counter, line_id_counter
            )

        if apply_cutlines:
            grid.set_obstacles_clones(db.get_obstacles())
        else:
            grid.set_obstacles_2d(db.get_obstacles())

        grid.set_boundary_conditions_2d(
            db.get_boundary_conditions_2d(),
            quadtree,
            node_id_counter,
            line_id_counter,
        )

        # This marks the nodes intersecting the dem_average geometries.
        dem_average_areas = db.get_dem_average_areas()
        grid.set_dem_averaged_cells(dem_average_areas)

    connection_nodes = db.get_connection_nodes()
    if grid_settings.use_1d_flow and len(connection_nodes) > 0:
        progress_callback(0.8, "Constructing 1D computational grid...")
        locations = db.get_cross_section_locations()
        pipes = db.get_pipes()
        weirs = db.get_weirs()
        culverts = db.get_culverts()
        orifices = db.get_orifices()
        cross_section_definitions = db.get_cross_section_definitions()
        (
            cross_section_definitions_unique,
            mapping,
        ) = cross_section_definitions.get_unique()
        map_cross_section_definition(
            [locations, orifices, pipes, culverts, weirs], mapping
        )

        cn_grid = Grid.from_connection_nodes(
            connection_nodes=connection_nodes,
            node_id_counter=node_id_counter,
        )
        connection_node_first_id = cn_grid.nodes.id[0] if len(cn_grid.nodes) > 0 else 0
        grid += cn_grid

        channels = db.get_channels()

        potential_breaches = db.get_potential_breaches()
        breach_points = potential_breaches.project_on_channels(channels).merge()

        channel_grid = Grid.from_linear_objects(
            connection_nodes=connection_nodes,
            objects=channels,
            fixed_nodes=breach_points,
            cell_tree=grid.cell_tree if grid_settings.use_2d else None,
            global_dist_calc_points=grid_settings.dist_calc_points,
            embedded_cutoff_threshold=grid_settings.embedded_cutoff_threshold,
            node_id_counter=node_id_counter,
            embedded_node_id_counter=embedded_node_id_counter,
            line_id_counter=line_id_counter,
            connection_node_offset=connection_node_first_id,
        )

        locations.apply_to_lines(channel_grid.lines, channels)
        windshieldings = db.get_windshieldings()
        windshieldings.apply_to_lines(channel_grid.lines)
        grid += channel_grid

        grid += Grid.from_linear_objects(
            connection_nodes=connection_nodes,
            objects=pipes,
            fixed_nodes=None,
            cell_tree=grid.cell_tree if grid_settings.use_2d else None,
            global_dist_calc_points=grid_settings.dist_calc_points,
            embedded_cutoff_threshold=grid_settings.embedded_cutoff_threshold,
            node_id_counter=node_id_counter,
            embedded_node_id_counter=embedded_node_id_counter,
            line_id_counter=line_id_counter,
            connection_node_offset=connection_node_first_id,
        )

        grid += Grid.from_linear_objects(
            connection_nodes=connection_nodes,
            objects=culverts,
            fixed_nodes=None,
            cell_tree=grid.cell_tree if grid_settings.use_2d else None,
            global_dist_calc_points=grid_settings.dist_calc_points,
            embedded_cutoff_threshold=grid_settings.embedded_cutoff_threshold,
            node_id_counter=node_id_counter,
            embedded_node_id_counter=embedded_node_id_counter,
            line_id_counter=line_id_counter,
            connection_node_offset=connection_node_first_id,
        )

        grid += grid.from_structures(
            connection_nodes=connection_nodes,
            weirs=weirs,
            orifices=orifices,
            line_id_counter=line_id_counter,
            connection_node_offset=connection_node_first_id,
        )
        grid.set_cross_sections(cross_section_definitions_unique)

        grid.set_calculation_types()
        grid.set_bottom_levels()
        grid.set_breach_ids(breach_points)
        grid.set_initial_waterlevels(
            connection_nodes=connection_nodes,
            channels=channels,
            pipes=pipes,
            culverts=culverts,
        )
        grid.set_boundary_conditions_1d(db.get_boundary_conditions_1d())
        grid.set_pumps(db.get_pumps())

    if grid.nodes.has_1d and grid.nodes.has_2d:
        progress_callback(0.9, "Connecting 1D and 2D domains...")
        grid.embed_nodes(embedded_node_id_counter)

        grid.add_1d2d_lines(
            exchange_lines=db.get_exchange_lines(),
            connection_nodes=connection_nodes,
            channels=channels,
            pipes=pipes,
            locations=locations,
            culverts=culverts,
            potential_breaches=potential_breaches,
            line_id_counter=line_id_counter,
            node_open_water_detection=grid_settings.node_open_water_detection,
        )
        if grid.meta.has_groundwater:
            channels.apply_has_groundwater_exchange(
                grid.nodes, grid.lines, grid.nodes_embedded
            )
            pipes.apply_has_groundwater_exchange(
                grid.nodes, grid.lines, grid.nodes_embedded
            )
            grid.add_1d2d_groundwater_lines(line_id_counter, connection_nodes)

    if grid_settings.use_0d_inflow in (
        InflowType.IMPERVIOUS_SURFACE.value,
        InflowType.SURFACE.value,
    ):
        # process zero-d information, either using the 'v2_surfaces'
        # or 'v2_impervious_surfaces' table.
        progress_callback(0.95, "Processing 0D domain...")
        surfaces: Surfaces = db.get_surfaces()
        grid.add_0d(surfaces)

    grid.finalize()
    return grid


def make_gridadmin(
    sqlite_path: Path,
    dem_path: Path,
    out_path: Optional[Path] = None,
    meta: dict = None,
    progress_callback: Optional[Callable[[float, str], None]] = None,
    upgrade: bool = False,
    convert_to_geopackage: bool = False,
    apply_cutlines: bool = False,
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
        upgrade: whether to upgrade the sqlite (inplace) before processing
        apply_cutlines: whether to apply (obstacles as) cutlines and create clone cells.

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
        upgrade=upgrade,
        convert_to_geopackage=convert_to_geopackage,
        apply_cutlines=apply_cutlines,
    )

    progress_callback(0.99, "Writing gridadmin...")
    with Writer(out_path) as writer:
        return writer.write(grid)


# Legacy
make_grid = make_gridadmin
