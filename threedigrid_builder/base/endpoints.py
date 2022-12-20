from typing import List, Optional

import numpy as np

from threedigrid_builder.constants import ContentType

from .array import Array
from .lines import Lines
from .nodes import Nodes

__all__ = ["Endpoints"]


class Endpoint:
    id: int
    line_idx: int
    is_start: bool
    node_id: int


class Endpoints(Array[Endpoint]):
    """Each line has 2 endpoints.

    This class makes it easier to reason about all the objects on a node.
    """

    scalars = ("lines",)

    @classmethod
    def for_connection_nodes(
        cls, nodes: Nodes, lines: Lines, line_types: Optional[List[ContentType]] = None
    ) -> "Endpoints":
        """Return the endpoints for connection nodes

        Optionally filter by the content type of the lines.
        """
        if line_types is None:
            where = np.arange(len(lines), dtype=int)
        else:
            where = np.where(np.isin(lines.content_type, line_types))[0]
        n = len(where)
        result = cls(
            lines=lines,
            id=range(n * 2),
            line_idx=np.repeat(where, 2),
            node_id=lines.line[where].ravel(),
            is_start=np.tile([True, False], n),
        )
        result = result[
            np.isin(
                result.node_id,
                nodes.id[nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES],
            )
        ]
        result.reorder_by("node_id")
        return result

    @property
    def is_end(self):
        return ~self.is_start

    @property
    def invert_level(self):
        return np.where(
            self.is_start, self.invert_level_start_point, self.invert_level_end_point
        )

    def __getattr__(self, name):
        return getattr(self.lines, name)[self.line_idx]

    def reduce_per_node(self, ufunc, values):
        """Compute the per-node ufunc reduction of 'values'

        For computing the maximum, use 'np.fmax'. For computing the sum, use 'np.add'.

        Returns node ids and corresponding maximum values.
        """
        if len(values) != len(self):
            raise ValueError("values must have the same length as self")
        if len(values) == 0:
            return np.empty(0, dtype=int), np.empty(0, dtype=float)
        diff = np.diff(self.node_id)
        if np.any(diff < 0):
            raise ValueError("node_id must be sorted")
        indices = np.concatenate([[0], np.where(diff > 0)[0] + 1])
        result = ufunc.reduceat(values, indices)
        return self.node_id[indices], result

    def nanmin_per_node(self, values):
        return self.reduce_per_node(np.fmin, values)
