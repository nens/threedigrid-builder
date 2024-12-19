from dataclasses import fields

import numpy as np
import shapely

from threedigrid_builder.base import (
    is_int_enum,
    is_tuple_type,
    Lines,
    OutputInterface,
    search,
    unpack_optional_type,
)
from threedigrid_builder.base.surfaces import SurfaceMaps, Surfaces
from threedigrid_builder.constants import ContentType, LineType, NodeType
from threedigrid_builder.grid import (
    Grid,
    GridMeta,
    Obstacles,
    PotentialBreaches,
    QuadtreeStats,
)
from threedigrid_builder.grid.cross_section_definitions import CrossSections

try:
    import h5py
except ImportError:
    h5py = None

__all__ = ["GridAdminOut"]


NODE_TYPES_1D = (NodeType.NODE_1D_NO_STORAGE, NodeType.NODE_1D_STORAGE)

LINE_TYPES_1D = (
    LineType.LINE_1D_EMBEDDED,
    LineType.LINE_1D_ISOLATED,
    LineType.LINE_1D_CONNECTED,
    LineType.LINE_1D_LONG_CRESTED,
    LineType.LINE_1D_SHORT_CRESTED,
    LineType.LINE_1D_DOUBLE_CONNECTED,
)

LINE_TYPES_1D2D = (
    LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED,
    LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER,
    LineType.LINE_1D2D_DOUBLE_CONNECTED_CLOSED,
    LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER,
)

LINE_TYPES_1D2D_GW = (LineType.LINE_1D2D_GROUNDWATER,)

LINE_TYPES_2D = (
    LineType.LINE_2D_U,
    LineType.LINE_2D_V,
    LineType.LINE_2D,
    LineType.LINE_2D_OBSTACLE,
    LineType.LINE_2D_OBSTACLE_U,
    LineType.LINE_2D_OBSTACLE_V,
)

# Some default settings for h5py datasets
HDF5_SETTINGS = {
    "compression": "gzip",  # more compatible than lzf and better compression
    "compression_opts": 1,  # the fastest, a 9 gives only 1% better compression
    "shuffle": True,  # helps another 3% and has almost no overhead
}


def field_to_h5(group, name, dtype, val, mode="attrs"):
    """Write a dataclass field to H5py attributes or datasets.

    None values are converted automatically to -9999 (int an IntEnum), NaN (float)
    and "" (str). bool values cannot be None. Unknown datatypes are skipped silently.
    """
    # backwards compat: epsg code is saved as string
    if name == "epsg_code":
        dtype = str
        val = str(val)
    # convert Optional[<type>] to <type>
    dtype = unpack_optional_type(dtype)
    # handle Tuple
    if is_tuple_type(dtype):
        if dtype.__args__[-1] is Ellipsis:
            shape = (len(val),)
        else:
            shape = (len(dtype.__args__),)
        dtype = dtype.__args__[0]
    else:
        shape = ()
    if dtype is bool:
        pass
    elif dtype is int or is_int_enum(dtype):
        dtype = np.int32
        if val is None:
            val = np.full(shape, -9999, dtype=dtype)
    elif dtype is float:
        dtype = np.float64
        if val is None:
            val = np.full(shape, np.nan, dtype=dtype)
    elif dtype is str:
        # Future: Write variable-length UTF-8
        # dtype = h5py.special_dtype(vlen=str)
        # But for calccore, we now need fixed-length ASCII
        if shape != ():
            raise RuntimeError("Cannot write string lists to HDF5.")
        if val is None:
            val = b""
        dtype = f"S{max(len(val), 1)}"
        if not isinstance(val, bytes):
            val = val.encode("utf-8", errors="ignore")
    else:
        return
    if mode == "attrs":
        group.attrs.create(name, data=val, shape=shape, dtype=dtype)
    elif mode == "datasets":
        group.create_dataset(name, data=val, shape=shape, dtype=dtype)
    else:
        raise ValueError(f"Unknown h5 write mode '{mode}'")


def dataclass_to_h5(group, datacls, mode="attrs"):
    """Write a dataclass to H5py attributes or datasets.

    None values are converted automatically to -9999 (int an IntEnum), NaN (float)
    and "" (str). bool values cannot be None. Unknown datatypes are skipped silently.
    """
    for field in fields(datacls):
        field_to_h5(group, field.name, field.type, getattr(datacls, field.name), mode)


def increase(arr):
    """Increase arr by one where arr is not -9999"""
    return np.add(arr, 1 * (arr != -9999), dtype=arr.dtype)


def to_bytes_array(arr, length):
    """Convert a list or array of strings to a numpy bytes array, utf-8 encoded."""
    result = np.zeros(len(arr), dtype=f"S{length}")
    for i, x in enumerate(arr):
        if x is None:
            continue
        elif isinstance(x, str):
            result[i] = x.encode("utf-8", errors="ignore")
        elif isinstance(x, bytes):
            result[i] = x
        else:
            raise TypeError(f"Unexpected type for bytes conversion, got '{type(x)}'")
    return result


class GridAdminOut(OutputInterface):
    def __init__(self, path):
        if not self.available():
            raise ImportError("Cannot write to HDF5 if h5py is not available.")
        super().__init__(path)

    @staticmethod
    def available():
        return h5py is not None

    def __enter__(self):
        self._file = h5py.File(self.path, mode="w")
        return self

    def __exit__(self, *args, **kwargs):
        self._file.close()

    def write(self, grid: Grid):
        self.write_meta(grid.meta)
        self.write_grid_counts(grid.nodes, grid.lines)
        self.write_quadtree(grid.quadtree_stats)
        self.write_nodes(grid.nodes)
        self.write_nodes_embedded(grid.nodes_embedded)
        self.write_quarters(grid.quarters)
        self.write_lines(grid.lines, grid.cross_sections)
        self.write_pumps(grid.pumps)
        self.write_cross_sections(grid.cross_sections)
        self.write_obstacles(grid.obstacles)
        self.write_breaches(grid.breaches)
        if grid.meta.has_0d:
            self.write_surfaces(grid.surfaces, grid.surface_maps)

    def write_meta(self, meta: GridMeta):
        """Write the metadata to the gridadmin file.

        Metadata is stored in attrs and settings are stored in dedicated groups.

        Args:
            meta (GridMeta)
        """
        dataclass_to_h5(self._file, meta, "attrs")

        grid_settings = self._file.create_group("grid_settings")
        dataclass_to_h5(grid_settings, meta.grid_settings, mode="datasets")

        tables_settings = self._file.create_group("tables_settings")
        dataclass_to_h5(tables_settings, meta.tables_settings, mode="datasets")

    def write_grid_counts(self, nodes, lines: Lines):
        """Write the "meta" group in the gridadmin file.

        Args:
            nodes (Nodes)
            lines (Lines)

        Raises:
            ValueError if it exists already.
        """
        group = self._file.create_group("meta")

        # the number of nodes in several categories
        for dataset_name, node_types in [
            ("n1dtot", NODE_TYPES_1D),
            ("n2dtot", (NodeType.NODE_2D_OPEN_WATER)),
            ("ngrtot", (NodeType.NODE_2D_GROUNDWATER)),
            ("n1dobc", (NodeType.NODE_1D_BOUNDARIES)),
            ("n2dobc", (NodeType.NODE_2D_BOUNDARIES)),
            ("ngr2bc", (NodeType.NODE_2D_GROUNDWATER_BOUNDARIES)),
        ]:
            count = np.count_nonzero(np.isin(nodes.node_type, node_types))
            group.create_dataset(dataset_name, data=count, dtype="i4")

        # the number of lines in several categories
        # note: remove the 1D boundary lines, these are counted via n1dobc
        masked_line_kcu = lines.kcu[lines.is_1d_boundary != 1]
        for dataset_name, kcu_values in [
            ("liutot", (LineType.LINE_2D_U, LineType.LINE_2D_OBSTACLE_U)),
            ("livtot", (LineType.LINE_2D_V, LineType.LINE_2D_OBSTACLE_V)),
            ("l1dtot", LINE_TYPES_1D),
            ("l2dtot", LINE_TYPES_2D),
            ("lgrtot", LineType.LINE_2D_VERTICAL),
            ("infl1d", LINE_TYPES_1D2D),
            ("ingrw1d", LINE_TYPES_1D2D_GW),
        ]:
            count = np.count_nonzero(np.isin(masked_line_kcu, kcu_values))
            group.create_dataset(dataset_name, data=count, dtype="i4")

        # Groundwater lines from open water lines
        for dataset_name, kcu_values in [
            ("lgutot", (LineType.LINE_2D_U, LineType.LINE_2D_OBSTACLE_U)),
            ("lgvtot", (LineType.LINE_2D_V, LineType.LINE_2D_OBSTACLE_V)),
        ]:
            if self._file.attrs["has_groundwater_flow"]:
                count = np.count_nonzero(np.isin(masked_line_kcu, kcu_values))
            else:
                count = 0
            group.create_dataset(dataset_name, data=count, dtype="i4")

        # the number of unique boundaries (only 2D)
        for dataset_name, node_type in [
            ("nob2ds", NodeType.NODE_2D_BOUNDARIES),
            ("nob2dg", NodeType.NODE_2D_GROUNDWATER_BOUNDARIES),
        ]:
            count = len(np.unique(nodes.boundary_id[nodes.node_type == node_type]))
            group.create_dataset(dataset_name, data=count, dtype="i4")

    def write_quadtree(self, quadtree_statistics):
        """Write the "grid_coordinate_attributes" group in the gridadmin file.

        Args:
            quadtree (QuadTree)

        Raises:
            ValueError if it exists already.
        """
        if quadtree_statistics is None:
            # Write a dummy for pure 1D models (Issue 217)
            quadtree_statistics = QuadtreeStats(
                lgrmin=0,
                kmax=1,
                mmax=[1],
                nmax=[1],
                dx=[-10000000.0],
                dxp=-1000000.0,
                x0p=0.0,
                y0p=0.0,
            )

        group = self._file.create_group("grid_coordinate_attributes")
        dataclass_to_h5(group, quadtree_statistics, mode="datasets")

    def write_nodes_embedded(self, nodes_embedded):
        """Write the "nodes_embedded" group in the gridadmin file

        See write_nodes.
        """
        return self.write_nodes(nodes_embedded, group_name="nodes_embedded")

    def write_nodes(self, nodes, group_name="nodes"):
        """Write the "nodes" group in the gridadmin file

        Raises a ValueError if it exists already.

        Notes:
            Some datasets were 64-bit integers, but now they are saved as 32-bit integers.
            The following datasets were added: code.
            For floats (double) we use NaN instead of -9999.0 to denote empty.

        Args:
            nodes (Nodes)
        """
        group = self._file.create_group(group_name)
        shape = (len(nodes),)

        # Some convenient masks:
        is_cn = nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES
        is_2d = np.isin(
            nodes.node_type, [NodeType.NODE_2D_OPEN_WATER, NodeType.NODE_2D_GROUNDWATER]
        )

        # Datasets that match directly to a nodes attribute:
        self.write_dataset(group, "id", nodes.id + 1)
        self.write_dataset(group, "code", to_bytes_array(nodes.code, 32))
        self.write_dataset(group, "node_type", nodes.node_type)
        self.write_dataset(group, "calculation_type", nodes.calculation_type)
        self.write_dataset(group, "coordinates", nodes.coordinates.T)
        self.write_dataset(group, "cell_coords", nodes.bounds.T)
        self.write_dataset(group, "nodk", nodes.nodk)
        self.write_dataset(group, "nodm", nodes.nodm)
        self.write_dataset(group, "nodn", nodes.nodn)
        self.write_dataset(group, "storage_area", nodes.storage_area)
        self.write_dataset(group, "dmax", nodes.dmax)
        self.write_dataset(group, "s1d", nodes.s1d)
        self.write_dataset(group, "embedded_in", nodes.embedded_in + 1)
        self.write_dataset(group, "boundary_id", nodes.boundary_id)
        self.write_dataset(group, "boundary_type", nodes.boundary_type)
        self.write_dataset(group, "has_dem_averaged", nodes.has_dem_averaged)
        self.write_dataset(group, "initial_waterlevel", nodes.initial_waterlevel)

        # content pk is only set for connection nodes, otherwise 0
        self.write_dataset(group, "content_pk", np.where(is_cn, nodes.content_pk, 0))
        # is_manhole is 1 for manholes, otherwise -9999
        self.write_dataset(
            group,
            "is_manhole",
            np.where(nodes.is_manhole, 1, -9999).astype("i4"),
        )
        self.write_dataset(group, "dimp", nodes.dimp)
        self.write_dataset(
            group, "bottom_level", np.where(nodes.is_manhole, nodes.dmax, np.nan)
        )
        self.write_dataset(group, "z_coordinate", nodes.dmax)

        # 2D stuff that is derived from 'bounds'
        pixel_width = np.full(len(nodes), -9999, dtype="i4")
        if is_2d.any():
            pixel_width[is_2d] = (
                nodes.pixel_coords[is_2d, 2] - nodes.pixel_coords[is_2d, 0]
            )

        self.write_dataset(group, "x_coordinate", nodes.coordinates[:, 0])
        self.write_dataset(group, "y_coordinate", nodes.coordinates[:, 1])
        self.write_dataset(group, "pixel_coords", nodes.pixel_coords.T)
        self.write_dataset(group, "pixel_width", pixel_width)
        self.write_dataset(
            group, "display_name", to_bytes_array(nodes.display_name, 64)
        )
        self.write_dataset(group, "zoom_category", nodes.zoom_category)
        # Set manhole_id to match nodes.id when there is a manhole
        self.write_dataset(
            group, "manhole_id", np.where(nodes.is_manhole, nodes.id + 1, -9999)
        )
        self.write_dataset(group, "manhole_indicator", nodes.manhole_indicator)
        self.write_dataset(group, "shape", to_bytes_array(nodes.shape, 4))
        self.write_dataset(group, "drain_level", nodes.drain_level)
        self.write_dataset(group, "surface_level", nodes.surface_level)
        self.write_dataset(group, "width", nodes.width)

        # filled in threedi-tables:
        self.write_dataset(group, "sumax", np.full(shape, np.nan, dtype=np.float64))

    def write_lines(self, lines: Lines, cross_sections):
        """Write the "lines" group in the gridadmin file

        Raises a ValueError if it exists already.

        Notes:
            Some datasets were 64-bit integers, but now they are saved as 32-bit integers.
            The following datasets were added: ds1d, dpumax, flod, flou, cross1,
            cross2, cross_weight
            For floats (double) we use NaN instead of -9999.0 to denote empty.

        Args:
            lines (Lines)
        """
        group = self._file.create_group("lines")
        shape = (len(lines),)
        is_channel = lines.content_type == ContentType.TYPE_V2_CHANNEL

        l2d = np.isin(lines.kcu, (LineType.LINE_2D_U, LineType.LINE_2D_V))
        l2d_obstacle = np.isin(
            lines.kcu, (LineType.LINE_2D_OBSTACLE_U, LineType.LINE_2D_OBSTACLE_V)
        )
        # Datasets that match directly to a lines attribute:
        self.write_dataset(group, "id", lines.id + 1)
        self.write_dataset(group, "code", to_bytes_array(lines.code, 32))
        self.write_dataset(
            group, "display_name", to_bytes_array(lines.display_name, 64)
        )

        lines.kcu[l2d] = LineType.LINE_2D
        lines.kcu[l2d_obstacle] = LineType.LINE_2D_OBSTACLE
        self.write_dataset(group, "kcu", lines.kcu)
        calculation_type = np.full(shape, -9999, dtype="i4")
        calculation_type[is_channel] = lines.kcu[is_channel] + 100
        self.write_dataset(group, "calculation_type", calculation_type)
        self.write_dataset(group, "line", lines.line.T + 1)
        self.write_dataset(group, "cross_pix_coords", lines.cross_pix_coords.T)
        self.write_dataset(group, "s1d", lines.s1d)
        self.write_dataset(group, "ds1d_half", lines.ds1d_half)
        self.write_dataset(group, "ds1d", lines.ds1d)
        self.write_dataset(group, "lik", lines.lik)
        self.write_dataset(group, "lim", lines.lim)
        self.write_dataset(group, "lin", lines.lin)

        content_type = np.full(len(lines), b"", dtype="S10")
        for ind in np.where(lines.content_type != -9999)[0]:
            content_type_name = ContentType(lines.content_type[ind]).name
            content_type[ind] = content_type_name.lstrip("TYPE_").lower()[:10]

        self.write_dataset(group, "content_type", content_type)
        self.write_dataset(group, "content_pk", lines.content_pk)
        self.write_dataset(group, "dpumax", lines.dpumax)
        self.write_dataset(
            group, "invert_level_start_point", lines.invert_level_start_point
        )
        self.write_dataset(
            group, "invert_level_end_point", lines.invert_level_end_point
        )
        self.write_dataset(group, "flod", lines.flod)
        self.write_dataset(group, "flou", lines.flou)

        cross1 = lines.cross_id1.copy()
        cross1[cross1 != -9999] = np.digitize(
            cross1[cross1 != -9999], cross_sections.content_pk, right=True
        )
        cross2 = lines.cross_id2.copy()
        cross2[cross2 != -9999] = np.digitize(
            cross2[cross2 != -9999], cross_sections.content_pk, right=True
        )
        self.write_dataset(group, "cross1", increase(cross1))
        self.write_dataset(group, "cross2", increase(cross2))
        self.write_dataset(group, "frict_type1", lines.frict_type1)
        self.write_dataset(group, "frict_type2", lines.frict_type2)
        self.write_dataset(group, "frict_value1", lines.frict_value1)
        self.write_dataset(group, "frict_value2", lines.frict_value2)
        self.write_dataset(group, "veg_coef1", lines.veg_coef1)
        self.write_dataset(group, "veg_coef2", lines.veg_coef2)
        self.write_dataset(group, "veg_height1", lines.veg_height1)
        self.write_dataset(group, "veg_height2", lines.veg_height2)
        self.write_dataset(group, "cross_weight", lines.cross_weight)
        self.write_dataset(
            group, "hydraulic_resistance_in", lines.hydraulic_resistance_in
        )
        self.write_dataset(
            group, "hydraulic_resistance_out", lines.hydraulic_resistance_out
        )
        self.write_dataset(group, "line_coords", lines.line_coords.T)
        self.write_dataset(
            group,
            "discharge_coefficient",
            lines.discharge_coefficient_positive,
        )
        self.write_dataset(
            group,
            "discharge_coefficient_negative",
            lines.discharge_coefficient_negative,
        )
        self.write_dataset(
            group,
            "discharge_coefficient_positive",
            lines.discharge_coefficient_positive,
        )

        # Transform an array of linestrings to list of coordinate arrays (x,x,y,y)
        self.write_line_geometry_dataset(
            group, "line_geometries", lines.line_geometries
        )
        self.write_dataset(group, "zoom_category", lines.zoom_category)
        self.write_dataset(
            group, "connection_node_end_pk", lines.connection_node_start_id
        )
        self.write_dataset(
            group, "connection_node_start_pk", lines.connection_node_end_id
        )
        self.write_dataset(group, "crest_level", lines.crest_level)
        self.write_dataset(group, "crest_type", lines.crest_type)

        # Cross section on pipes, culverts, weirs, and orifices
        pipe_culvert = np.logical_or.reduce(
            (
                lines.content_type == ContentType.TYPE_V2_PIPE,
                lines.content_type == ContentType.TYPE_V2_CULVERT,
                lines.content_type == ContentType.TYPE_V2_WEIR,
                lines.content_type == ContentType.TYPE_V2_ORIFICE,
            )
        )
        # Data placeholders
        cross_section_width = np.full(len(lines), np.nan, dtype=np.float64)
        cross_section_height = np.full(len(lines), np.nan, dtype=np.float64)
        cross_section_shape = np.full(len(lines), -9999, dtype="i4")

        # Cross section definition ids, for pipe and culvert cross_id1 == cross_id2
        cross_ids = lines.cross_id1[pipe_culvert]

        # Put cross section data into placeholders at pipe and culvert indexes
        # Get cross section width, height, shape corresponding to cross section definition id
        np.put(
            cross_section_width,
            lines.id[pipe_culvert],
            cross_sections.width_1d[
                search(
                    cross_sections.content_pk, cross_ids, mask=None, assume_ordered=True
                )
            ],
        )
        np.put(
            cross_section_height,
            lines.id[pipe_culvert],
            cross_sections.height_1d[
                search(
                    cross_sections.content_pk, cross_ids, mask=None, assume_ordered=True
                )
            ],
        )
        np.put(
            cross_section_shape,
            lines.id[pipe_culvert],
            cross_sections.shape[
                search(
                    cross_sections.content_pk, cross_ids, mask=None, assume_ordered=True
                )
            ],
        )
        self.write_dataset(group, "cross_section_width", cross_section_width)
        self.write_dataset(group, "cross_section_height", cross_section_height)
        self.write_dataset(group, "cross_section_shape", cross_section_shape)
        self.write_dataset(group, "dist_calc_points", lines.dist_calc_points)
        self.write_dataset(group, "material", lines.material)
        self.write_dataset(group, "sewerage", lines.sewerage)
        self.write_dataset(group, "sewerage_type", lines.sewerage_type)
        has_friction_data = np.isin(
            lines.content_pk,
            [
                ContentType.TYPE_V2_PIPE,
                ContentType.TYPE_V2_CULVERT,
                ContentType.TYPE_V2_WEIR,
                ContentType.TYPE_V2_ORIFICE,
            ],
        )
        self.write_dataset(
            group,
            "friction_type",
            np.where(has_friction_data, lines.frict_type1, -9999),
        )
        self.write_dataset(
            group,
            "friction_value",
            np.where(has_friction_data, lines.frict_value1, np.nan),
        )
        self.write_dataset(group, "windshieldings", lines.windshieldings.T)

    def write_quarters(self, quarters, group_name="quarters"):
        """Write the "quarters" group in the gridadmin file

        Raises a ValueError if it exists already.

        Notes:
            Some datasets were 64-bit integers, but now they are saved as 32-bit integers.

        Args:
            quarters (Quarters)
        """
        if quarters is None:
            return
        group = self._file.create_group(group_name)
        mask = quarters.id != -9999
        quarters.id[mask] = quarters.id[mask] + 1
        quarters.id[~mask] = 0
        self.write_dataset(group, "id", quarters.id)

        mask = quarters.line != -9999
        quarters.line[mask] = quarters.line[mask] + 1
        quarters.line[~mask] = 0
        self.write_dataset(group, "line", quarters.line.T)

        mask = quarters.neighbour_node != -9999
        quarters.neighbour_node[mask] = quarters.neighbour_node[mask] + 1
        quarters.neighbour_node[~mask] = 0
        self.write_dataset(group, "neighbour_node", quarters.neighbour_node.T)

    def write_pumps(self, pumps):
        group = self._file.create_group("pumps")

        # Datasets that match directly to a lines attribute:
        self.write_dataset(group, "bottom_level", pumps.bottom_level)
        self.write_dataset(group, "code", to_bytes_array(pumps.code, 32))
        self.write_dataset(group, "content_pk", pumps.content_pk)
        self.write_dataset(group, "id", pumps.id + 1)
        self.write_dataset(group, "node1_id", increase(pumps.line[:, 0]))
        self.write_dataset(group, "node2_id", increase(pumps.line[:, 1]))
        self.write_dataset(group, "node_coordinates", pumps.line_coords.T)
        self.write_dataset(group, "capacity", pumps.capacity)
        self.write_dataset(
            group, "connection_node_end_pk", pumps.connection_node_end_id
        )
        self.write_dataset(
            group, "connection_node_start_pk", pumps.connection_node_start_id
        )
        self.write_dataset(group, "coordinates", pumps.coordinates.T)
        self.write_dataset(group, "lower_stop_level", pumps.lower_stop_level)
        self.write_dataset(group, "start_level", pumps.start_level)
        self.write_dataset(group, "type", pumps.type_)
        self.write_dataset(group, "upper_stop_level", pumps.upper_stop_level)
        self.write_dataset(
            group, "display_name", to_bytes_array(pumps.display_name, 64)
        )
        self.write_dataset(group, "zoom_category", pumps.zoom_category)

    def write_surfaces(self, surfaces: Surfaces, surface_maps: SurfaceMaps):
        # For now use 0 instead of NaN as fill value for surfaces
        # The calcore does not support NaN (for surfaces) yet
        default_fill_value = 0

        group = self._file.create_group("surface")
        self.write_dataset(group, "id", surfaces.id, fill=default_fill_value)
        self.write_dataset(group, "code", to_bytes_array(surfaces.code, 100))
        self.write_dataset(
            group, "display_name", to_bytes_array(surfaces.display_name, 250)
        )
        self.write_dataset(group, "function", to_bytes_array(surfaces.function, 64))

        self.write_dataset(
            group,
            "area",
            surfaces.area,
            fill=default_fill_value,
            fill_nan=default_fill_value,
        )
        self.write_dataset(
            group, "centroid_x", surfaces.centroid_x, fill=default_fill_value
        )
        self.write_dataset(
            group, "centroid_y", surfaces.centroid_y, fill=default_fill_value
        )
        self.write_dataset(
            group,
            "dry_weather_flow",
            surfaces.dry_weather_flow,
            fill=default_fill_value,
            fill_nan=default_fill_value,
        )
        self.write_dataset(
            group,
            "nr_of_inhabitants",
            surfaces.nr_of_inhabitants,
            fill=default_fill_value,
            fill_nan=default_fill_value,
        )
        self.write_dataset(group, "infiltration_flag", surfaces.infiltration_flag)
        self.write_dataset(
            group,
            "outflow_delay",
            surfaces.outflow_delay,
            fill=default_fill_value,
            fill_nan=default_fill_value,
        )
        self.write_dataset(
            group,
            "storage_limit",
            surfaces.storage_limit,
            fill=default_fill_value,
            fill_nan=default_fill_value,
        )

        if surfaces.surface_class is not None and np.any(
            surfaces.surface_class != None  # NOQA
        ):  # noqa
            # Impervious surfaces
            self.write_dataset(
                group, "surface_class", to_bytes_array(surfaces.surface_class, 128)
            )
            self.write_dataset(
                group,
                "surface_inclination",
                to_bytes_array(surfaces.surface_inclination, 64),
            )
            self.write_dataset(
                group,
                "surface_sub_class",
                to_bytes_array(surfaces.surface_sub_class, 128),
            )

        # Surface params
        self.write_dataset(group, "fb", surfaces.fb, fill=default_fill_value)
        self.write_dataset(group, "fe", surfaces.fe, fill=default_fill_value)
        self.write_dataset(group, "ka", surfaces.ka, fill=default_fill_value)
        self.write_dataset(group, "kh", surfaces.kh, fill=default_fill_value)

        # Surface maps (surface to connectionodes)
        self.write_dataset(group, "fac", surface_maps.fac, fill=default_fill_value)
        self.write_dataset(group, "nxc", surface_maps.nxc, fill=default_fill_value)
        self.write_dataset(group, "nyc", surface_maps.nyc, fill=default_fill_value)
        self.write_dataset(group, "pk", surface_maps.pk, fill=default_fill_value)

        # Note: +1 for Fortran 1-based indexing
        self.write_dataset(
            group, "cci", increase(surface_maps.cci), fill=default_fill_value
        )
        self.write_dataset(
            group, "imp", increase(surface_maps.imp), fill=default_fill_value
        )

    def write_cross_sections(self, cross_sections: CrossSections):
        if len(cross_sections) == 0:
            return
        group = self._file.create_group("cross_sections")

        # Datasets that match directly to a lines attribute:
        self.write_dataset(group, "id", cross_sections.id + 1)
        self.write_dataset(group, "code", to_bytes_array(cross_sections.code, 32))
        self.write_dataset(group, "shape", cross_sections.shape)
        self.write_dataset(group, "content_pk", cross_sections.content_pk)
        self.write_dataset(group, "width_1d", cross_sections.width_1d)
        self.write_dataset(group, "offset", cross_sections.offset)
        self.write_dataset(group, "count", cross_sections.count)
        self.write_dataset(group, "tables", cross_sections.tables.T, insert_dummy=False)
        if cross_sections.tables_yz is not None:
            self.write_dataset(group, "offset_yz", cross_sections.offset_yz)
            self.write_dataset(group, "count_yz", cross_sections.count_yz)
            self.write_dataset(
                group, "tables_yz", cross_sections.tables_yz.T, insert_dummy=False
            )

    def write_obstacles(self, obstacles: Obstacles):
        """For backwards compat, the group is named 'levees'"""
        if obstacles is None:
            return
        group = self._file.create_group("levees")

        self.write_dataset(group, "id", obstacles.id, insert_dummy=False)
        self.write_dataset(
            group, "crest_level", obstacles.crest_level, insert_dummy=False
        )
        self.write_dataset(
            group,
            "max_breach_depth",
            np.full(len(obstacles), fill_value=np.nan),
            insert_dummy=False,
        )
        self.write_line_geometry_dataset(
            group, "coords", obstacles.the_geom, insert_dummy=False
        )

    def write_breaches(self, breaches: PotentialBreaches):
        if breaches is None:
            return
        group = self._file.create_group("breaches")

        self.write_dataset(group, "id", breaches.id + 1)
        self.write_dataset(group, "levl", breaches.line_id + 1)
        self.write_dataset(group, "content_pk", breaches.content_pk)
        self.write_dataset(group, "levbr", breaches.maximum_breach_depth)
        self.write_dataset(group, "levmat", breaches.levee_material)
        self.write_dataset(
            group, "coordinates", shapely.get_coordinates(breaches.the_geom).T
        )
        self.write_dataset(group, "code", to_bytes_array(breaches.code, 32))
        self.write_dataset(
            group, "display_name", to_bytes_array(breaches.display_name, 64)
        )

    def write_dataset(
        self, group, name, values, fill=None, insert_dummy=True, fill_nan=None
    ):
        """Create the correct size dataset for writing to gridadmin.h5 and
        filling the extra indices with correct fillvalues.

        Args:
            group (hdf5 group): group to write to in hdf5 file.
            name (str): Name of dataset
            values (array): Values of dataset to write.
            fill: fillvalue with same dtype as values.
            fill_nane: (optional) fillvalue for nan values
        """

        if fill_nan is not None:
            if np.issubdtype(values.dtype, np.floating):
                values = values.copy()
                values[np.isnan(values)] = fill_nan

        if fill is None:
            if np.issubdtype(values.dtype, np.floating):
                fill = np.nan
            elif np.issubdtype(values.dtype, bytes):
                fill = b""
            else:
                fill = -9999
        assert insert_dummy in (True, False)

        if values.ndim == 1:
            ds = group.create_dataset(
                name,
                (values.shape[0] + int(insert_dummy),),
                dtype=values.dtype,
                fillvalue=fill,
                **HDF5_SETTINGS,
            )
            ds[int(insert_dummy) :] = values
            if insert_dummy and name == "id":  # set the ID of the dummy element
                ds[0] = 0
        elif values.ndim == 2:
            ds = group.create_dataset(
                name,
                (values.shape[0], values.shape[1] + int(insert_dummy)),
                dtype=values.dtype,
                fillvalue=fill,
                **HDF5_SETTINGS,
            )
            ds[:, int(insert_dummy) :] = values
        else:
            ValueError("Too many dimensions for values.")

    def write_line_geometry_dataset(self, group, name, data, insert_dummy=True):
        # Transform an array of linestrings to list of coordinate arrays (x,x,y,y)
        line_geometries = [shapely.get_coordinates(x).T.ravel() for x in data]
        if insert_dummy:
            line_geometries.insert(0, np.array([np.nan, np.nan]))
        # The dataset has a special "variable length" dtype
        try:
            vlen_dtype = h5py.vlen_dtype(np.dtype(float))
        except AttributeError:  # Pre h5py 2.10
            vlen_dtype = h5py.special_dtype(vlen=np.dtype(float))

        # insert line geometry data preserving its original type
        geometry_data = np.empty(len(line_geometries), dtype=object)
        geometry_data[:] = line_geometries

        group.create_dataset(
            name,
            data=geometry_data,
            dtype=vlen_dtype,
            **HDF5_SETTINGS,
        )
