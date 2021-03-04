from threedigrid_builder.base import OutputInterface
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import NodeType

import numpy as np


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

    def write_nodes(self, nodes, pixel_size, **kwargs):
        """Write the "nodes" group in the gridadmin file

        Raises a ValueError if it exists already.

        Notes:
            The following datasets (node attributes) where previously 64-bit, but now
            32-bit integer: calculation_type, is_manhole, manhole_indicator,
            pixel_coords, pixel_width, zoom_category.
            For floats (double precision) we use NaN instead of -9999.0 to denote empty.

        Args:
            nodes (Nodes)
            pixel_size (float): the size of a pixel in projected units (mostly meters)
        """
        group = self._file.create_group("nodes")

        # Some convenient masks:
        is_mh = nodes.content_type == ContentType.TYPE_V2_MANHOLE
        is_cn = is_mh | (nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES)
        is_2d = nodes.node_type == NodeType.NODE_2D_OPEN_WATER

        # Datasets that match directly to a nodes attribute:
        group.create_dataset("id", data=nodes.id)
        group.create_dataset("node_type", data=nodes.node_type)
        group.create_dataset("calculation_type", data=nodes.calculation_type)
        group.create_dataset("coordinates", data=nodes.coordinates.T)
        group.create_dataset("cell_coords", data=nodes.bounds.T)
        group.create_dataset("storage_area", data=nodes.storage_area)

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
        # missing in gridadmin: code
        group.create_dataset(
            "display_name", data=np.full(len(nodes), -9999, dtype="S64")
        )
        group.create_dataset(
            "zoom_category", data=np.full(len(nodes), -9999, dtype="i4")
        )
        group.create_dataset("initial_waterlevel", len(nodes), dtype=float)
        # (manhole specific:)
        group.create_dataset(
            "manhole_indicator", data=np.full(len(nodes), -9999, dtype="i4")
        )
        group.create_dataset("shape", data=np.full(len(nodes), b"-999", dtype="S4"))
        group.create_dataset("drain_level", len(nodes), dtype=float)
        group.create_dataset("surface_level", len(nodes), dtype=float)
        group.create_dataset("width", len(nodes), dtype=float)  # but no length?

        # unknown
        group.create_dataset("seq_id", data=nodes.id)
        group.create_dataset("sumax", len(nodes), dtype=float)

    def write_lines(self, lines, **kwargs):
        pass
