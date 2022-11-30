from typing import Tuple

import numpy as np
import pygeos

from threedigrid_builder.base import Array, Nodes, search
from threedigrid_builder.constants import CalculationType, ContentType
from threedigrid_builder.exceptions import SchematisationError

__all__ = ["ExchangeLines"]


class ExchangeLine:
    id: int
    the_geom: pygeos.Geometry
    channel_id: int


class ExchangeLines(Array[ExchangeLine]):
    def split_primary_secondary(self) -> Tuple[ExchangeLine, ExchangeLine]:
        """Split the exchange lines into primary and secondary ones"""
        exc1_id, exc1_idx, counts = np.unique(
            self.channel_id, return_index=True, return_counts=True
        )
        if np.any(counts > 2):
            raise SchematisationError(
                f"Channels {exc1_id.tolist()} have too many exchange lines."
            )
        exc2_idx = np.delete(np.arange(len(self)), exc1_idx)
        return self[exc1_idx], self[exc2_idx]

    def find_for_nodes(self, content_pk, content_type):
        """Return an exchange line id (or -9999 if not present) or each node"""
        result = np.full(len(content_pk), fill_value=-9999, dtype=self.id.dtype)
        if len(self) == 0:
            return result
        is_channel = content_type == ContentType.TYPE_V2_CHANNEL
        idx = search(
            self.channel_id,
            content_pk[is_channel],
            assume_ordered=False,
            check_exists=False,
        )
        ids = self.id[idx]
        ids[self.channel_id[idx] != content_pk[is_channel]] = -9999
        result[is_channel] = ids
        return result

    def get_line_mappings(self, nodes: Nodes):
        """Return a node index and an exchange line id per 1D2D line.

        Returns 2 arrays of shape (n_lines, ):
        - an int array with a node index
        - an int array with a exchange line id (-9999 if there is none)
        """
        primary, secondary = self.split_primary_secondary()
        # Initialize the arrays
        is_single = nodes.calculation_type == CalculationType.CONNECTED
        is_double = nodes.calculation_type == CalculationType.DOUBLE_CONNECTED
        node_idx_1 = np.where(is_single | is_double)[0]
        node_idx_2 = np.where(is_double)[0]

        # for the channels, insert exchange line idx
        exc_id_1 = primary.find_for_nodes(
            nodes.content_pk[node_idx_1], nodes.content_type[node_idx_1]
        )
        exc_id_2 = secondary.find_for_nodes(
            nodes.content_pk[node_idx_2], nodes.content_type[node_idx_2]
        )

        # Merge the arrays
        node_idx = np.concatenate([node_idx_1, node_idx_2])
        exc_id = np.concatenate([exc_id_1, exc_id_2])

        # Sort them by node index
        sorter = np.argsort(node_idx)
        return node_idx[sorter], exc_id[sorter]

    def get_2d_side_points(self, nodes: Nodes):
        """Return a 2D side Point for each line

        Returns 3 arrays of shape (n_lines, ):
        - an int array with a node index
        - an int array with a exchange line id (-9999 if there is none)
        - a Geometry array with points
        """
        # Create an array of node indices & corresponding exchange line indices
        node_idx, exc_id = self.get_line_mappings(nodes)
        has_exc = exc_id != -9999

        # Collect the 1D sides of the 1D2D line
        geometries = pygeos.points(nodes.coordinates[node_idx])

        # Change these into the 2D sides through the exchange line
        # if there is no exchange line, 2D side equals 1D side.
        exchange_line_geoms = self.the_geom[self.id_to_index(exc_id[has_exc])]
        pygeos.prepare(exchange_line_geoms)
        geometries[has_exc] = pygeos.get_point(
            pygeos.shortest_line(exchange_line_geoms, geometries[has_exc]), 0
        )

        return node_idx, exc_id, geometries

    def get_lines_node_idx(self, nodes: Nodes, cell_tree: pygeos.STRtree):
        """Return a 1D and 2D node index per line

        Returns 3 arrays of shape (n_lines, ):
        - an int array with a 2D node (cell) index
        - an int array with a 1D node index
        - an int array with a exchange line id (-9999 if there is none)
        """
        node_idx, exc_id, side_2d = self.get_2d_side_points(nodes)
        # The query_bulk returns 2 1D arrays: one with indices into the supplied node
        # geometries and one with indices into the tree of cells.
        idx = cell_tree.query_bulk(side_2d)
        # Address edge cases of multiple 1D-2D lines per node: just take the one
        _, first_unique_index = np.unique(idx[0], return_index=True)
        idx = idx[:, first_unique_index]
        n_lines = idx.shape[1]
        # Error if there is a node without a 1D-2D line
        if n_lines != node_idx.shape[0]:
            # The code in this if clause is only for pretty error formatting.
            out_of_bounds = np.delete(node_idx, idx[0])
            out_of_bounds_exc = np.delete(exc_id, idx[0])
            object_pk_list = _get_pk_content_type(
                nodes, out_of_bounds[out_of_bounds_exc == -9999]
            )
            msg = ""
            if len(object_pk_list) > 0:
                msg += (
                    f"The following object(s) have a connected calculation type but are "
                    f"(partially) outside of the 2D calculation cells: "
                    f"{', '.join(object_pk_list)}. "
                )
            if (out_of_bounds_exc != -9999).any():
                mask = out_of_bounds_exc != -9999
                exc_ids_list = out_of_bounds_exc[mask].tolist()
                msg += (
                    f"The following exchange line(s) are outside of a 2D calculation "
                    f"cell: {str(exc_ids_list)}."
                )
            raise SchematisationError(msg)

        return idx[1], node_idx[idx[0]], exc_id


HAS_1D2D_CONTENT_TYPES = (
    ContentType.TYPE_V2_CONNECTION_NODES,
    ContentType.TYPE_V2_CHANNEL,
    ContentType.TYPE_V2_PIPE,
    ContentType.TYPE_V2_CULVERT,
)


def _get_pk_content_type(nodes, out_of_bounds):
    """Determine the object pk and object type of calculation points for pretty
    printing.

    Args:
      nodes (Nodes):
      out_of_bounds (ndarray): Array with node idx of nodes that are out of bounds.

    Returns:
      List of content types (their pretty name) and pks of objects with nodes
      out of bounds.
    """
    if len(out_of_bounds) == 0:
        return []
    types = nodes.content_type[out_of_bounds]
    pks = nodes.content_pk[out_of_bounds]
    pretty_names = ("connection nodes", "channels", "pipes", "culverts")
    return [
        f"{pretty_name} {sorted(set(pks[types == content_type]))}"
        for content_type, pretty_name in zip(HAS_1D2D_CONTENT_TYPES, pretty_names)
        if np.any(types == content_type)
    ]
