from dataclasses import fields
from threedigrid_builder.base import is_int_enum
from threedigrid_builder.base import is_tuple_type
from threedigrid_builder.base import OutputInterface
from threedigrid_builder.base import unpack_optional_type
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.grid import GridMeta
from threedigrid_builder.grid.cross_section_definitions import CrossSections

import numpy as np
import pygeos


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
    LineType.LINE_1D2D_POSSIBLE_BREACH,
    LineType.LINE_1D2D_ACTIVE_BREACH,
)

LINE_TYPES_1D2D_GW = (
    LineType.LINE_1D2D_GROUNDWATER_OPEN_WATER,
    LineType.LINE_1D2D_GROUNDWATER_SEWER,
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
            val = val.encode()
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

    def write_grid_counts(self, nodes, lines):
        """Write the "meta" group in the gridadmin file.

        Args:
            nodes (Nodes)
            lines (Lines)

        Raises:
            ValueError if it exists already.
        """
        group = self._file.create_group("meta")
        NODATA_YET = 0

        n2dtot = np.count_nonzero(nodes.node_type == NodeType.NODE_2D_OPEN_WATER)
        n2dobc = np.count_nonzero(nodes.node_type == NodeType.NODE_2D_BOUNDARIES)
        n1dobc = np.count_nonzero(nodes.node_type == NodeType.NODE_1D_BOUNDARIES)
        group.create_dataset("n2dtot", data=n2dtot, dtype="i4")
        group.create_dataset("n2dobc", data=n2dobc, dtype="i4")
        group.create_dataset("ngr2bc", data=NODATA_YET, dtype="i4")

        n1dtot = np.count_nonzero(np.isin(nodes.node_type, NODE_TYPES_1D))
        group.create_dataset("n1dtot", data=n1dtot, dtype="i4")
        group.create_dataset("n1dobc", data=n1dobc, dtype="i4")
        liutot = np.count_nonzero(lines.kcu == LineType.LINE_2D_U) + np.count_nonzero(
            lines.kcu == LineType.LINE_2D_OBSTACLE_U
        )
        group.create_dataset("liutot", data=liutot, dtype="i4")
        livtot = np.count_nonzero(lines.kcu == LineType.LINE_2D_V) + np.count_nonzero(
            lines.kcu == LineType.LINE_2D_OBSTACLE_V
        )
        group.create_dataset("livtot", data=livtot, dtype="i4")

        group.create_dataset("lgutot", data=NODATA_YET, dtype="i4")
        group.create_dataset("lgvtot", data=NODATA_YET, dtype="i4")

        l1dtot = np.count_nonzero(np.isin(lines.kcu, LINE_TYPES_1D))
        group.create_dataset("l1dtot", data=l1dtot, dtype="i4")

        infl1d = np.count_nonzero(np.isin(lines.kcu, LINE_TYPES_1D2D))
        group.create_dataset("infl1d", data=infl1d, dtype="i4")

        ingrw1d = np.count_nonzero(np.isin(lines.kcu, LINE_TYPES_1D2D_GW))
        group.create_dataset("ingrw1d", data=ingrw1d, dtype="i4")
        group.create_dataset("jap1d", data=NODATA_YET, dtype="i4")

    def write_quadtree(self, quadtree_statistics):
        """Write the "grid_coordinate_attributes" group in the gridadmin file.

        Args:
            quadtree (QuadTree)

        Raises:
            ValueError if it exists already.
        """

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
        is_2d = nodes.node_type == NodeType.NODE_2D_OPEN_WATER

        # Datasets that match directly to a nodes attribute:
        self.write_dataset(group, "id", nodes.id + 1)
        self.write_dataset(group, "code", nodes.code.astype("S32"), fill=b"")
        self.write_dataset(group, "display_name", nodes.code.astype("S64"), fill=b"")
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
        self.write_dataset(group, "initial_waterlevel", nodes.initial_waterlevel)

        # content pk is only set for connection nodes, otherwise 0
        self.write_dataset(group, "content_pk", np.where(is_cn, nodes.content_pk, 0))
        # is_manhole is 1 for manholes, otherwise -9999
        self.write_dataset(
            group,
            "is_manhole",
            np.where(nodes.manhole_id != -9999, 1, -9999).astype("i4"),
        )

        # unclear what is the difference: seems bottom_level is for 1D, and z_coordinate
        # for 2D, but some nodes have both sets (then they are equal)
        # missing in gridadmin: dimp (bottom level groundwater)
        self.write_dataset(group, "bottom_level", nodes.dmax)
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

        # can be collected from SQLite, but empty for now:
        self.write_dataset(
            group, "zoom_category", np.full(len(nodes), -9999, dtype="i4")
        )
        # (manhole specific:)
        self.write_dataset(
            group, "manhole_indicator", np.full(len(nodes), -9999, dtype="i4")
        )
        self.write_dataset(
            group, "shape", np.full(len(nodes), b"-999", dtype="S4"), fill=b""
        )
        self.write_dataset(
            group, "drain_level", np.full(shape, -9999, dtype=np.float64)
        )
        self.write_dataset(
            group, "surface_level", np.full(shape, -9999, dtype=np.float64)
        )
        self.write_dataset(group, "width", np.full(shape, -9999, dtype=np.float64))

        # unknown
        self.write_dataset(group, "sumax", np.full(shape, -9999, dtype=np.float64))

    def write_lines(self, lines):
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
        fill_int = np.full(shape, -9999, dtype="i4")
        fill_float = np.full(shape, -9999, dtype=float)
        is_channel = lines.content_type == ContentType.TYPE_V2_CHANNEL

        l2d = np.isin(lines.kcu, (LineType.LINE_2D_U, LineType.LINE_2D_V))
        l2d_obstacle = np.isin(
            lines.kcu, (LineType.LINE_2D_OBSTACLE_U, LineType.LINE_2D_OBSTACLE_V)
        )
        # Datasets that match directly to a lines attribute:
        self.write_dataset(group, "id", lines.id + 1)
        self.write_dataset(group, "code", lines.code.astype("S32"), fill=b"")
        self.write_dataset(
            group, "display_name", lines.display_name.astype("S64"), fill=b""
        )

        lines.kcu[l2d] = LineType.LINE_2D
        lines.kcu[l2d_obstacle] = LineType.LINE_2D_OBSTACLE
        lines.kcu[lines.kcu == LineType.LINE_1D_BOUNDARY] = LineType.LINE_1D_ISOLATED
        self.write_dataset(group, "kcu", lines.kcu)
        calculation_type = fill_int.copy()
        calculation_type[is_channel] = lines.kcu[is_channel] + 100
        self.write_dataset(group, "calculation_type", calculation_type)
        self.write_dataset(group, "line", lines.line.T + 1)
        self.write_dataset(group, "cross_pix_coords", lines.cross_pix_coords.T)
        self.write_dataset(group, "s1d", lines.s1d)
        self.write_dataset(group, "ds1d", lines.ds1d)
        self.write_dataset(group, "lik", lines.lik)
        self.write_dataset(group, "lim", lines.lim)
        self.write_dataset(group, "lin", lines.lin)

        content_type = np.full(len(lines), b"", dtype="S10")
        for ind in np.where(lines.content_type != -9999)[0]:
            content_type_name = ContentType(lines.content_type[ind]).name
            content_type[ind] = content_type_name.lstrip("TYPE_").lower()[:10]

        self.write_dataset(group, "content_type", content_type, fill=b"")
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
        self.write_dataset(group, "cross1", increase(lines.cross1))
        self.write_dataset(group, "cross2", increase(lines.cross2))
        self.write_dataset(group, "frict_type1", lines.frict_type1)
        self.write_dataset(group, "frict_type2", lines.frict_type2)
        self.write_dataset(group, "frict_value1", lines.frict_value1)
        self.write_dataset(group, "frict_value2", lines.frict_value2)
        self.write_dataset(group, "cross_weight", lines.cross_weight)
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
        line_geometries = [
            pygeos.get_coordinates(x).T.ravel() for x in lines.line_geometries
        ]
        line_geometries.insert(0, np.array([-9999.0, -9999.0]))
        # The dataset has a special "variable length" dtype. This one is special write method.
        try:
            vlen_dtype = h5py.vlen_dtype(np.dtype(float))
        except AttributeError:  # Pre h5py 2.10
            vlen_dtype = h5py.special_dtype(vlen=np.dtype(float))
        group.create_dataset(
            "line_geometries",
            data=np.array(line_geometries, dtype=object),
            dtype=vlen_dtype,
            **HDF5_SETTINGS,
        )

        # can be collected from SQLite, but empty for now:
        self.write_dataset(group, "connection_node_end_pk", fill_int)
        self.write_dataset(group, "connection_node_start_pk", fill_int)
        self.write_dataset(group, "crest_level", fill_float)
        self.write_dataset(group, "crest_type", fill_int)
        self.write_dataset(group, "cross_section_height", fill_int)
        self.write_dataset(group, "cross_section_shape", fill_int)
        self.write_dataset(group, "cross_section_width", fill_float)
        self.write_dataset(group, "dist_calc_points", fill_float)
        self.write_dataset(group, "friction_type", fill_int)
        self.write_dataset(group, "friction_value", fill_float)
        self.write_dataset(group, "material", fill_int)
        self.write_dataset(group, "sewerage", fill_int)
        self.write_dataset(group, "sewerage_type", fill_int)
        self.write_dataset(group, "zoom_category", fill_int)

    def write_pumps(self, pumps):
        group = self._file.create_group("pumps")

        # Datasets that match directly to a lines attribute:
        self.write_dataset(group, "bottom_level", pumps.bottom_level)
        self.write_dataset(group, "code", pumps.code.astype("S32"), fill=b"")
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

        # can be collected from SQLite, but empty for now:
        self.write_dataset(
            group, "display_name", np.full(len(pumps), b"", dtype="S64"), fill=b""
        )
        self.write_dataset(
            group, "zoom_category", np.full(len(pumps), -9999, dtype="i4")
        )

    def write_cross_sections(self, cross_sections: CrossSections):
        group = self._file.create_group("cross_sections")

        # Datasets that match directly to a lines attribute:
        self.write_dataset(group, "id", cross_sections.id + 1)
        self.write_dataset(group, "code", cross_sections.code.astype("S32"), fill=b"")
        self.write_dataset(group, "shape", cross_sections.shape)
        self.write_dataset(group, "content_pk", cross_sections.content_pk)
        self.write_dataset(group, "width_1d", cross_sections.width_1d)
        self.write_dataset(group, "offset", cross_sections.offset)
        self.write_dataset(group, "count", cross_sections.count)

        # do not use self.write_dataset as we don't want a dummy element
        group.create_dataset("tables", data=cross_sections.tables.T, **HDF5_SETTINGS)

    def write_dataset(self, group, name, values, fill=None):
        """Create the correct size dataset for writing to gridadmin.h5 and
        filling the extra indices with correct fillvalues.

        Args:
            group (hdf5 group): group to write to in hdf5 file.
            name (str): Name of dataset
            values (array): Values of dataset to write.
            fill: fillvalue with same dtype as values.
        """
        if fill is None:
            if np.issubdtype(values.dtype, np.floating):
                fill = np.nan
            else:
                fill = -9999

        if values.ndim == 1:
            ds = group.create_dataset(
                name,
                (values.shape[0] + 1,),
                dtype=values.dtype,
                fillvalue=fill,
                **HDF5_SETTINGS,
            )
            ds[1:] = values
            if name == "id":  # set the ID of the dummy element
                ds[0] = 0
        elif values.ndim == 2:
            ds = group.create_dataset(
                name,
                (values.shape[0], values.shape[1] + 1),
                dtype=values.dtype,
                fillvalue=fill,
                **HDF5_SETTINGS,
            )
            ds[:, 1:] = values
        else:
            ValueError("Too many dimensions for values.")
