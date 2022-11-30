from typing import Tuple

import numpy as np
import pygeos

from threedigrid_builder.base import Array, search
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

    def get_line_mappings(self, nodes):
        """Return a node index and an exchange line id per 1D2D line.

        Returns 3 arrays of shape (n_lines, ):
        - an int array with a node index
        - an int array with a exchange line index (-9999 if there is none)
        - a bool array with True where a line is double connected
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
