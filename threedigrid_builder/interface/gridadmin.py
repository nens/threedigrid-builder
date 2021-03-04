import h5py
import numpy as np


class GridAdmin:
    def __init__(self, filename):
        # create the file if it does not exist
        self._file = h5py.File(filename, 'a')

    def write_nodes(self, nodes, pixel_size, origin):
        """Write the "nodes" group in the gridadmin file

        Raises a ValueError if it exists already.

        Args:
            nodes (Nodes)
            pixel_size (float): the size of a pixel in projected units (mostly meters)
            origin (tuple of 2 floats): the coordinate of pixel (0, 0)
        """
        group = self._file.create_group("nodes")

        # converted the following datasets from 64 to 32-bit integer (i8 to i4):
        # - calculation_type
        # - is_manhole
        # - manhole_indicator
        # - pixel_coords
        # - pixel_width
        # - zoom_category
        # replaced -9999 in floats to NaN

        # id and seq_id are the same??
        group.create_dataset("id", data=nodes.id)
        group.create_dataset("seq_id", data=nodes.id)

        group.create_dataset("node_type", data=nodes.node_type)
        group.create_dataset("calculation_type", data=nodes.calculation_type)
        # missing in gridadmin: content_type
        group.create_dataset("content_pk", data=nodes.content_pk)
        group.create_dataset("coordinates", data=nodes.coordinates.T)
        group.create_dataset("cell_coords", data=nodes.bounds.T)
        group.create_dataset("bottom_level", data=nodes.dmax)  # only for 1D

        # only for 2D? some 1D nodes also have this, but then it is equal to bottom_level
        group.create_dataset("z_coordinate", data=nodes.dmax)

        # missing in gridadmin: dimp (bottom level groundwater)
        group.create_dataset("storage_area", data=nodes.storage_area)

        # can be collected from SQLite, but empty for now:
        group.create_dataset("display_name", data=np.full(len(nodes), -9999, dtype="S64"))
        group.create_dataset("zoom_category", data=np.full(len(nodes), -9999, dtype="i4"))
        group.create_dataset("is_manhole", data=np.full(len(nodes), -9999, dtype="i4"))
        group.create_dataset("manhole_indicator", data=np.full(len(nodes), -9999, dtype="i4"))
        group.create_dataset("shape", data=np.full(len(nodes), b"-999", dtype="S4"))

        # not fully clear 2D stuff
        x = np.full(len(nodes), np.nan, dtype=float)
        y = np.full(len(nodes), np.nan, dtype=float)
        is_2d = np.isfinite(nodes["cell_coords"][0])
        x[is_2d] = (nodes.bounds[is_2d, 0] + nodes.bounds[is_2d, 2]) / 2
        y[is_2d] = (nodes.bounds[is_2d, 1] + nodes.bounds[is_2d, 3]) / 2
        pixel_coords = np.full((4, len(nodes)), -9999, dtype="i4")
        pixel_coords[:, is_2d] = np.round(
            (nodes.bounds[is_2d] - (tuple(origin) * 2)).T / pixel_size
        )
        pixel_width = np.full(len(nodes), -9999, dtype="i4")
        pixel_width[is_2d] = pixel_coords[2] - pixel_coords[0]
        group.create_dataset("x_coordinate", data=x)
        group.create_dataset("y_coordinate", data=y)
        group.create_dataset("z_coordinate", data=nodes.dmax)  # ??
        group.create_dataset("pixel_coords", data=pixel_coords)
        group.create_dataset("pixel_width", data=pixel_width)

        # unknown
        group.create_dataset("drain_level", len(nodes), dtype=float)
        group.create_dataset("initial_waterlevel", len(nodes), dtype=float)
        group.create_dataset("sumax", len(nodes), dtype=float
        group.create_dataset("surface_level", len(nodes), dtype=float)
        group.create_dataset("width", len(nodes), dtype=float)
