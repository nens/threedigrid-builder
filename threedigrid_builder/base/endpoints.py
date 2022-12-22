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
        cls,
        nodes: Nodes,
        lines: Lines,
        line_mask=None,
        node_mask=None,
    ) -> "Endpoints":
        """Return the endpoints for connection nodes

        Optionally filter by the content type of the lines.
        """
        if line_mask is None:
            where = np.arange(len(lines), dtype=int)
        else:
            where = np.where(line_mask)[0]
        n = len(where)
        result = cls(
            lines=lines,
            id=range(n * 2),
            line_idx=np.repeat(where, 2),
            node_id=lines.line[where].ravel(),
            is_start=np.tile([True, False], n),
        )
        is_cn = nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES
        node_mask = is_cn if node_mask is None else is_cn & node_mask
        return result[np.isin(result.node_id, nodes.id[node_mask])]

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

    def reduce_per_node(self, reduceat, values) -> "NodeValues":
        """Compute the per-node reduction of 'values'

        The 'reduceat' function is normally a numpy ufunc.reduceat. This function
        is expected to reduce 'values' on index ranges contained in 'indices'. The length
        of 'indices' must equal the length of the returned value.

        For instance:

         - self.node_id = [1, 1, 5, 5, 5]
         - values = [1., 3., 5., 7., 9.]
         - gives: indices = [0, 2]
         - and reduceat is expected to act on indices [0, 2) and [2, <end>)

        Returns node ids and corresponding maximum values.
        """
        if len(values) != len(self):
            raise ValueError("values must have the same length as self")
        if len(values) == 0:
            return NodeValues(id=[])
        diff = np.diff(self.node_id)
        if np.any(diff < 0):
            raise ValueError("endpoints must be ordered by node_id")
        indices = np.concatenate([[0], np.where(diff > 0)[0] + 1])
        result = reduceat(values, indices)
        return NodeValues(id=self.node_id[indices], value=result)

    def nanmin_per_node(self, values) -> "NodeValues":
        return self.reduce_per_node(np.fmin.reduceat, values)

    def first_per_node(self, values) -> "NodeValues":
        return self.reduce_per_node(lambda x, y: x[y], values)


class NodeValue:
    id: int
    value: float


class NodeValues(Array[NodeValue]):
    pass
