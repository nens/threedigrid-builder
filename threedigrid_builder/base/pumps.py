from threedigrid_builder.base import array_of
from threedigrid_builder.constants import ContentType
from typing import Tuple

import numpy as np


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
    # zoom_category
    # display_name
    # sewerage
    # Generated fields:
    content_pk: int
    line: Tuple[int, int]
    line_coords: Tuple[float, float, float, float]
    coordinates: Tuple[float, float]
    bottom_level: float


@array_of(Pump)
class Pumps:
    def renumber(self):
        """Renumber the ids of the pumps with a sequential ID.

        The original ID is stored in content_pk.
        """
        self.content_pk[:] = self.id
        self.id[:] = np.arange(len(self))

    def set_node_data(self, nodes):
        """Set the node ids into self.line by mapping connection_node_start/end_id."""
        mask = nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES
        node_ids = nodes.id[mask]
        connection_node_ids = nodes.content_pk[mask]

        # sort by connection node id so that we can lookup
        sorter = np.argsort(connection_node_ids)
        node_ids = node_ids[sorter]
        connection_node_ids = connection_node_ids[sorter]

        self.line[:, 0] = node_ids[
            np.searchsorted(connection_node_ids, self.connection_node_start_id)
        ]
        has_end = self.connection_node_end_id >= 0
        self.line[has_end, 1] = node_ids[
            np.searchsorted(connection_node_ids, self.connection_node_end_id[has_end])
        ]

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
