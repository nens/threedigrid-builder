from threedigrid_builder.base import array_of, Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ManholeIndicator, ContentType, NodeType

import numpy as np
import pygeos
import itertools

__all__ = ["ConnectionNodes"]


class ConnectionNode:
    id: int
    the_geom: pygeos.Geometry
    code: str
    storage_area: float
    manhole_id: int
    calculation_type: CalculationType
    manhole_indicator: ManholeIndicator
    bottom_level: float
    drain_level: float
    surface_level: float
    manhole_shape: str  # enum with classes "00", "01", "02"
    manhole_width: float


@array_of(ConnectionNode)
class ConnectionNodes:
    def get_nodes(self, node_id_counter):
        """Convert connection nodes to a nodes instance

        Args:
            node_id_counter (iterable): an iterable yielding integers
            global_dist_calc_points (float): Default node interdistance.

        Returns:
            tuple of nodes (Nodes), segment_size (ndarray)
            nodes has data in the following columns:
            - id: 0-based counter generated here
            - coordinates
            - content_type: TYPE_V2_CHANNEL
            - content_pk: the id of the Channel from which this node originates
            - node_type: NODE_1D_STORAGE or NODE_1D_NO_STORAGE depending on storage_area
        """
        node_type = np.where(
            self.storage_area > 0,
            int(NodeType.NODE_1D_STORAGE),
            int(NodeType.NODE_1D_NO_STORAGE),
        )
        nodes = Nodes(
            id=itertools.islice(node_id_counter, len(self)),
            coordinates=pygeos.get_coordinates(self.the_geom),
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            content_pk=self.id,
            node_type=node_type,
        )
        return nodes
