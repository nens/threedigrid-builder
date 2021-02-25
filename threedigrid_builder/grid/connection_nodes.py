from .grid import Grid  # TODO this import should not be here
from threedigrid_builder.base import array_of
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import ContentType

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
            - nodes.content_pk: the user-supplied id
        """
        nodes = Nodes(
            id=range(len(self)),
            coordinates=pygeos.get_coordinates(self.the_geom),
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            content_pk=self.id,
        )
        lines = Lines(id=[])
        return Grid(nodes, lines)
