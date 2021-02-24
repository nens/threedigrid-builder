from threedigrid_builder.base import array_of
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import Grid
from threedigrid_builder.grid import Lines
from threedigrid_builder.grid import Nodes

import numpy as np
import pygeos


__all__ = ["ConnectionNodes"]


class ConnectionNode:
    id: int
    the_geom: pygeos.Geometry
    code: str
    storage_area: float


@array_of(ConnectionNode)
class ConnectionNodes:
    def get_grid(self):
        """Compute the grid (only nodes) for the connection nodes

        Args:
            self (ConnectionNodes):
                id and the_geom are used

        Returns:
            Grid with data in the following columns:
            - nodes.id: 0-based counter generated here
            - nodes.coordinates: node coordinates (from self.the_geom)
            - nodes.content_type: ContentType.TYPE_V2_CONNECTION_NODES
            - nodes.content_pk: connection node primary key (from self.id)
        """
        nodes = Nodes(
            id=np.arange(self.id.size),
            coordinates=pygeos.get_coordinates(self.the_geom),
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            content_pk=self.id,
        )
        lines = Lines(
            id=[]
        )
        return Grid(nodes, lines)
