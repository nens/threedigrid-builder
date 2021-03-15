from threedigrid_builder.base import OutputInterface
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.constants import LineType

import numpy as np
import pygeos


try:
    import h5py
except ImportError:
    h5py = None

__all__ = ["GridAdminOut"]


class GridAdminOut(OutputInterface):
    def __init__(self, path):
        if h5py is None:
            raise ImportError("Cannot write to HDF5 if h5py is not available.")
        super().__init__(path)

    def __enter__(self):
        self._file = h5py.File(self.path, mode="w")
        return self

    def __exit__(self, *args, **kwargs):
        self._file.close()

    def write_attrs(self, nodes, lines, epsg_code):
        """Write the "attrs" of the gridadmin file.

        Raises a ValueError if it exists already.
        Args:
            attrs of grid. (dict)
        """

        self._file.attrs.create("epsg_code", epsg_code)

        is_1d = (self.nodes.node_type == NodeType.NODE_1D_NO_STORAGE) + (
            self.nodes.node_type == NodeType.NODE_1D_STORAGE
        )
        if is_1d.any():
            extent_1d = np.array(
                [
                    np.amin(self.nodes.coordinates[is_1d, 0]),
                    np.amin(self.nodes.coordinates[is_1d, 1]),
                    np.amax(self.nodes.coordinates[is_1d, 0]),
                    np.amax(self.nodes.coordinates[is_1d, 1]),
                ]
            )
            self._file.attrs.create("extent_1d", extent_1d)
            self._file.attrs.create("has_1d", 1)
        else:
            self._file.attrs.create(
                "extent_1d", np.array([-9999.0, -9999.0, -9999.0, -9999.0])
            )
            self._file.attrs.create("has_1d", 0)

        is_2d = self.nodes.node_type == NodeType.NODE_2D_OPEN_WATER
        if is_2d.any:
            extent_2d = np.array(
                [
                    np.amin(self.nodes.coordinates[is_2d, 0]),
                    np.amin(self.nodes.coordinates[is_2d, 1]),
                    np.amax(self.nodes.coordinates[is_2d, 0]),
                    np.amax(self.nodes.coordinates[is_2d, 1]),
                ]
            )
            self._file.attrs.create("extent_2d", extent_2d)
            self._file.attrs.create("has_2d", 1)
        else:
            self._file.attrs.create(
                "extent_2d", np.array([-9999.0, -9999.0, -9999.0, -9999.0])
            )
            self._file.attrs.create("has_2d", 0)
        self._file.attrs.create("has_interception", 0)
        self._file.attrs.create("has_pumpstations", 0)
        self._file.attrs.create("has_simple_infiltration", 0)
        self._file.attrs.create("model_name", "...")
        self._file.attrs.create("model_slug", "...")
        self._file.attrs.create("revision_hash", "...")
        self._file.attrs.create("revision_nr", 0)
        self._file.attrs.create("threedigrid_builder_version", 0)

    def write_meta(self, meta):
        """Write the "meta" group in the gridadmin file.

        Raises a ValueError if it exists already.
        Args:
            meta attribute of grid. (dict)
        """
        group = self._file.create_group("meta")
        NODATA_YET = 0

        n2dtot = np.count_nonzero(self.nodes.node_type == NodeType.NODE_2D_OPEN_WATER)
        group.create_dataset("n2dtot", data=n2dtot, dtype="i4")
        group.create_dataset("n2dobc", data=NODATA_YET, dtype="i4")
        group.create_dataset("ngr2bc", data=NODATA_YET, dtype="i4")

        n1dtot = np.count_nonzero(
            self.nodes.node_type == NodeType.NODE_1D_NO_STORAGE
        ) + np.count_nonzero(self.nodes.node_type == NodeType.NODE_1D_STORAGE)
        group.create_dataset("n1dtot", data=n1dtot, dtype="i4")
        group.create_dataset("n1dobc", data=NODATA_YET, dtype="i4")

        liutot = np.count_nonzero(self.lines.line_type == LineType.LINE_2D_U)
        group.create_dataset("liutot", data=liutot, dtype="i4")
        livtot = np.count_nonzero(self.lines.line_type == LineType.LINE_2D_V)
        group.create_dataset("livtot", data=livtot, dtype="i4")

        group.create_dataset("lgutot", data=NODATA_YET, dtype="i4")
        group.create_dataset("lgvtot", data=NODATA_YET, dtype="i4")

        line_types_1d = [
            LineType.LINE_1D_EMBEDDED,
            LineType.LINE_1D_ISOLATED,
            LineType.LINE_1D_CONNECTED,
            LineType.LINE_1D_LONG_CRESTED,
            LineType.LINE_1D_SHORT_CRESTED,
            LineType.LINE_1D_DOUBLE_CONNECTED,
        ]
        l1dtot = sum(
            [np.count_nonzero(self.lines.line_type == x) for x in line_types_1d]
        )
        group.create_dataset("l1dtot", data=l1dtot, dtype="i4")

        line_types_1d2d = [
            LineType.LINE_1D2D_SINGLE_CONNECTED_WITH_STORAGE,
            LineType.LINE_1D2D_SINGLE_CONNECTED_WITHOUT_STORAGE,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_WITH_STORAGE,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_WITHOUT_STORAGE,
            LineType.LINE_1D2D_POSSIBLE_BREACH,
            LineType.LINE_1D2D_ACTIVE_BREACH,
        ]
        infl1d = sum(
            [np.count_nonzero(self.lines.line_type == x) for x in line_types_1d2d]
        )
        group.create_dataset("infl1d", data=infl1d, dtype="i4")

        line_types_1d2d_gw = [
            LineType.LINE_1D2D_GROUNDWATER_57,
            LineType.LINE_1D2D_GROUNDWATER_58,
        ]
        ingrw1d = sum(
            [np.count_nonzero(self.lines.line_type == x) for x in line_types_1d2d_gw]
        )
        group.create_dataset("ingrw1d", data=ingrw1d, dtype="i4")
        group.create_dataset("jap1d", data=NODATA_YET, dtype="i4")

    def write_quadtree(self, quadtree):
        """Write the "grid_coordinate_attributes" group in the gridadmin file.

        Raises a ValueError if it exists already.
        Args:
            quadtree (QuadTree)
        """

        group = self._file.create_group("grid_coordinate_attributes")
        group.create_dataset("lgrmin", data=quadtree.lgrmin, dtype="i4")
        group.create_dataset("dx", data=quadtree.dx)
        group.create_dataset("kmax", data=quadtree.kmax, dtype="i4")
        group.create_dataset("mmax", data=quadtree.mmax, dtype="i4")
        group.create_dataset("nmax", data=quadtree.nmax, dtype="i4")
        group.create_dataset("x0p", data=quadtree.bbox[0])
        group.create_dataset("y0p", data=quadtree.bbox[1])

    def write_nodes(self, nodes, pixel_size, **kwargs):
        """Write the "nodes" group in the gridadmin file

        Raises a ValueError if it exists already.

        Notes:
            Some datasets were 64-bit integers, but now they are saved as 32-bit integers.
            The following datasets were added: code.
            For floats (double) we use NaN instead of -9999.0 to denote empty.

        Args:
            nodes (Nodes)
            pixel_size (float): the size of a pixel in projected units (mostly meters)
        """
        group = self._file.create_group("nodes")
        shape = (len(nodes),)

        # Some convenient masks:
        is_mh = nodes.content_type == ContentType.TYPE_V2_MANHOLE
        is_cn = is_mh | (nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES)
        is_2d = nodes.node_type == NodeType.NODE_2D_OPEN_WATER

        # Datasets that match directly to a nodes attribute:
        group.create_dataset("id", data=nodes.id)
        group.create_dataset("code", data=nodes.code.astype("S32"))
        group.create_dataset("display_name", data=nodes.code.astype("S64"))
        group.create_dataset("node_type", data=nodes.node_type)
        group.create_dataset("calculation_type", data=nodes.calculation_type)
        group.create_dataset("coordinates", data=nodes.coordinates.T)
        group.create_dataset("cell_coords", data=nodes.bounds.T)
        group.create_dataset("nodk", data=nodes.nodk)
        group.create_dataset("nodm", data=nodes.nodm)
        group.create_dataset("nodn", data=nodes.nodn)
        group.create_dataset("storage_area", data=nodes.storage_area.astype("S32"))

        # content pk is only set for connection nodes, otherwise 0
        content_pk = np.full(len(nodes), 0, dtype="i4")
        content_pk[is_cn] = nodes.content_pk[is_cn]
        group.create_dataset("content_pk", data=content_pk)
        # is_manhole is 1 for manholes, otherwise -9999
        is_manhole = is_mh.astype("i4")
        is_manhole[~is_mh] = -9999
        group.create_dataset("is_manhole", data=is_manhole)

        # unclear what is the difference: seems bottom_level is for 1D, and z_coordinate
        # for 2D, but some nodes have both sets (then they are equal)
        # missing in gridadmin: dimp (bottom level groundwater)
        group.create_dataset("bottom_level", data=nodes.dmax)
        group.create_dataset("z_coordinate", data=nodes.dmax)

        # 2D stuff that is derived from 'bounds'
        x_coordinate = np.full(len(nodes), np.nan, dtype=float)
        y_coordinate = np.full(len(nodes), np.nan, dtype=float)
        pixel_coords = np.full((4, len(nodes)), -9999, dtype="i4")
        pixel_width = np.full(len(nodes), -9999, dtype="i4")
        if is_2d.any():
            x_coordinate[is_2d] = (nodes.bounds[is_2d, 0] + nodes.bounds[is_2d, 2]) / 2
            y_coordinate[is_2d] = (nodes.bounds[is_2d, 1] + nodes.bounds[is_2d, 3]) / 2
            origin = tuple(nodes.bounds[is_2d, :2].min(axis=0))
            pixel_coords[:, is_2d] = np.round(
                (nodes.bounds[is_2d] - (origin * 2)).T / pixel_size
            )
            pixel_width[is_2d] = pixel_coords[2] - pixel_coords[0]

        group.create_dataset("x_coordinate", data=x_coordinate)
        group.create_dataset("y_coordinate", data=y_coordinate)
        group.create_dataset("pixel_coords", data=pixel_coords)
        group.create_dataset("pixel_width", data=pixel_width)

        # can be collected from SQLite, but empty for now:
        group.create_dataset(
            "zoom_category", data=np.full(len(nodes), -9999, dtype="i4")
        )
        group.create_dataset("initial_waterlevel", shape, dtype=float)
        # (manhole specific:)
        group.create_dataset(
            "manhole_indicator", data=np.full(len(nodes), -9999, dtype="i4")
        )
        group.create_dataset("shape", data=np.full(len(nodes), b"-999", dtype="S4"))
        group.create_dataset("drain_level", shape, dtype=float)
        group.create_dataset("surface_level", shape, dtype=float)
        group.create_dataset("width", shape, dtype=float)  # but no length?

        # unknown
        group.create_dataset("sumax", shape, dtype=float)

    def write_lines(self, lines, **kwargs):
        """Write the "lines" group in the gridadmin file

        Raises a ValueError if it exists already.

        Notes:
            Some datasets were 64-bit integers, but now they are saved as 32-bit integers.
            The following datasets were added: ds1d, dpumax, flod, flou, cross1, cross2,
            cross_weight
            For floats (double) we use NaN instead of -9999.0 to denote empty.

        Args:
            lines (Lines)
        """
        group = self._file.create_group("lines")
        shape = (len(lines),)

        l2d = lines.line_type == LineType.LINE_2D_U
        l2d += lines.line_type == LineType.LINE_2D_V
        # Datasets that match directly to a lines attribute:
        group.create_dataset("id", data=lines.id)
        group.create_dataset("code", data=lines.code.astype("S32"))
        group.create_dataset("display_name", data=lines.display_name.astype("S64"))

        lines.line_type[l2d] = LineType.LINE_2D
        group.create_dataset("kcu", data=lines.line_type)
        group.create_dataset("calculation_type", data=lines.calculation_type)
        group.create_dataset("line", data=lines.line.T)
        group.create_dataset("ds1d", data=lines.ds1d)
        group.create_dataset("lik", data=lines.lik)
        group.create_dataset("lim", data=lines.lim)
        group.create_dataset("lin", data=lines.lin)

        content_type = np.full(len(lines), b"", dtype="S10")
        for ind in np.where(lines.content_type != -9999)[0]:
            content_type_name = ContentType(lines.content_type[ind]).name
            content_type[ind] = content_type_name.lstrip("TYPE_").lower()[:10]

        group.create_dataset("content_type", data=content_type)
        group.create_dataset("content_pk", data=lines.content_pk)
        group.create_dataset("dpumax", data=lines.dpumax)
        group.create_dataset("flod", data=lines.flod)
        group.create_dataset("flou", data=lines.flou)
        group.create_dataset("cross1", data=lines.cross1)
        group.create_dataset("cross2", data=lines.cross2)
        group.create_dataset("cross_weight", data=lines.cross_weight)
        group.create_dataset("line_coords", data=lines.line_coords.T)

        # Transform an array of linestrings to list of coordinate arrays (x,y,x,y,...)
        coords = pygeos.get_coordinates(lines.line_geometries)
        end = np.cumsum(pygeos.get_num_coordinates(lines.line_geometries))
        start = np.roll(end, 1)
        start[0] = 0
        line_geometries = [coords[a:b].ravel() for a, b in zip(start, end)]
        # The dataset has a special "variable length" dtype
        vlen_dtype = h5py.special_dtype(vlen=np.dtype(float))
        group.create_dataset("line_geometries", data=line_geometries, dtype=vlen_dtype)

        # can be collected from SQLite, but empty for now:
        group.create_dataset("connection_node_end_pk", shape, dtype="i4")
        group.create_dataset("connection_node_start_pk", shape, dtype="i4")
        group.create_dataset("crest_level", shape, dtype=float)
        group.create_dataset("crest_type", shape, dtype="i4")
        group.create_dataset("cross_section_height", shape, dtype="i4")
        group.create_dataset("cross_section_shape", shape, dtype="i4")
        group.create_dataset("cross_section_width", shape, dtype=float)
        group.create_dataset("discharge_coefficient", shape, dtype=float)
        group.create_dataset("discharge_coefficient_negative", shape, dtype=float)
        group.create_dataset("discharge_coefficient_positive", shape, dtype=float)
        group.create_dataset("dist_calc_points", shape, dtype=float)
        group.create_dataset("friction_type", shape, dtype="i4")
        group.create_dataset("friction_value", shape, dtype=float)
        group.create_dataset("invert_level_end_point", shape, dtype=float)
        group.create_dataset("invert_level_start_point", shape, dtype=float)
        group.create_dataset("material", shape, dtype="i4")
        group.create_dataset("sewerage", shape, dtype="i4")
        group.create_dataset("sewerage_type", shape, dtype="i4")
        group.create_dataset("zoom_category", shape, dtype="i4")
