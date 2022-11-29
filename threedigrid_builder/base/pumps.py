from typing import Tuple

import numpy as np

from threedigrid_builder.base import Array
from threedigrid_builder.constants import ContentType
from threedigrid_builder.exceptions import SchematisationError

__all__ = ["Pumps"]


class Pump:
    id: int
    code: str
    capacity: float
    connection_node_start_id: int
    connection_node_end_id: int  # optional!
    type_: int
    start_level: float
    lower_stop_level: float
    upper_stop_level: float
    zoom_category: int
    display_name: str
    # sewerage
    # Generated fields:
    content_pk: int
    line: Tuple[int, int]
    line_coords: Tuple[float, float, float, float]
    coordinates: Tuple[float, float]
    bottom_level: float


class Pumps(Array[Pump]):
    def renumber(self):
        """Renumber the ids of the pumps with a sequential ID.

        The original ID is stored in content_pk.
        """
        self.content_pk[:] = self.id
        self.id[:] = np.arange(len(self))

    def set_lines(self, nodes):
        """Set the node ids into self.line by mapping connection_node_start/end_id."""
        mask = nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES
        node_ids = nodes.id[mask]
        cn_ids = nodes.content_pk[mask]

        if (np.diff(cn_ids) < 0).any():
            # sort by connection node id so that we can lookup
            sorter = np.argsort(cn_ids)
            node_ids = node_ids[sorter]
            cn_ids = cn_ids[sorter]

        # check which connection node ids are specified
        has_start = self.connection_node_start_id >= 0
        has_end = self.connection_node_end_id >= 0
        if not has_start.all():
            raise SchematisationError(
                f"Pumps {self.content_pk[~has_start].tolist()} have no "
                f"connection_node_start_id"
            )

        # lookup the corresponding index into node_ids
        start_idx = np.searchsorted(cn_ids, self.connection_node_start_id)
        end_idx = np.searchsorted(cn_ids, self.connection_node_end_id[has_end])

        # check if the connection node ids indeed exist
        start_is_ok = (
            cn_ids.take(start_idx, mode="clip") == self.connection_node_start_id
        )
        end_is_ok = (
            cn_ids.take(end_idx, mode="clip") == self.connection_node_end_id[has_end]
        )
        if not (start_is_ok.all() and end_is_ok.all()):
            raise SchematisationError(
                f"Pumps {self.content_pk[~(start_is_ok & end_is_ok)].tolist()} refer "
                f"to nonexisting connection nodes."
            )
        self.line[:, 0] = node_ids[start_idx]
        self.line[has_end, 1] = node_ids[end_idx]

    def set_node_data(self, nodes):
        """Set bottom_level, line_coords, and coordinates from node ids in self.line."""
        has_end = self.line[:, 1] != -9999
        start_node_idx = nodes.id_to_index(self.line[:, 0])
        end_node_idx = nodes.id_to_index(self.line[has_end, 1])

        # copy over the bottom levels
        self.bottom_level[:] = nodes.dmax[start_node_idx]

        # set the line coordinates
        self.line_coords[:, :2] = nodes.coordinates[start_node_idx]
        self.line_coords[has_end, 2:] = nodes.coordinates[end_node_idx]

        # set the coordinates
        self.coordinates[:] = self.line_coords[:, :2]
        # take the centroid for pumps that have 2 coordinates
        self.coordinates[has_end] = (
            self.coordinates[has_end] + self.line_coords[has_end, 2:]
        ) / 2
