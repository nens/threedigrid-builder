import numpy as np

from threedigrid_builder.constants import ContentType

from .array import Array
from .lines import Lines
from .nodes import Nodes

__all__ = ["LineHalfs"]


class LineHalf:
    id: int
    line_idx: int
    is_start: bool
    node_id: int


class LineHalfs(Array[LineHalf]):
    """Each line has 2 line_halfs.

    This class makes it easier to reason about all the objects on a node.
    """

    scalars = ("lines",)

    @classmethod
    def from_nodes_lines(
        cls,
        nodes: Nodes,
        lines: Lines,
        line_mask=None,
        node_mask=None,
    ) -> "LineHalfs":
        """Return the line_halfs, optionally filtering with node / line masks."""
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
        if node_mask is not None:
            result = result[np.isin(result.node_id, nodes.id[node_mask])]
        return result

    @classmethod
    def for_connection_nodes(
        cls,
        nodes: Nodes,
        lines: Lines,
        line_mask=None,
        node_mask=None,
    ) -> "LineHalfs":
        """Return the line_halfs for connection nodes

        Optionally filter by the content type of the lines.
        """
        is_cn = nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES
        node_mask = is_cn if node_mask is None else is_cn & node_mask
        return cls.from_nodes_lines(nodes, lines, line_mask, node_mask)

    @property
    def is_end(self):
        return ~self.is_start

    @property
    def invert_level(self):
        return np.where(
            self.is_start, self.invert_level_start_point, self.invert_level_end_point
        )

    @property
    def ds1d(self):
        return np.where(
            self.is_start,
            self.lines.ds1d_half[self.line_idx],
            self.lines.ds1d[self.line_idx] - self.lines.ds1d_half[self.line_idx],
        )

    def __getattr__(self, name):
        try:
            return super(self.__class__, self).__getattr__(name)
        except AttributeError:
            return getattr(self.lines, name)[self.line_idx]

    def _get_reduce_indices(self):
        """Return indices into self containing the same node ids.

        Intended for usage in reduction functions (see reduce).

        For instance:

         - self.node_id = [1, 1, 5, 5, 5]
         - gives: indices = [0, 2]
        """
        if len(self) == 0:
            return np.empty((0,), dtype=int)
        diff = np.diff(self.node_id)
        if np.any(diff < 0):
            raise ValueError("line_halfs must be ordered by node_id")
        return np.concatenate([[0], np.where(diff > 0)[0] + 1])

    def get_reduce_id(self):
        """Return unique node ids.

        Intended for usage with reduction functions (see reduce).
        """
        return self.node_id[self._get_reduce_indices()]

    def reduce(self, reduceat, values):
        """Compute the per-node reduction of 'values'

        The 'reduceat' function is normally a numpy ufunc.reduceat. This function
        is expected to reduce 'values' on index ranges contained in 'indices'. The length
        of 'indices' must equal the length of the returned value.

        For instance:

         - self.node_id = [1, 1, 5, 5, 5]
         - values = [1., 3., 5., 7., 9.]
         - gives: indices = [0, 2]
         - and reduceat is expected to act on index ranges [0, 2) and [2, <end>)

        Returns node ids and corresponding maximum values.
        """
        values = np.asarray(values)
        if values.shape[0] != len(self):
            raise ValueError("values must have the same length as self")
        if len(values) == 0:
            return np.empty((0,) + values.shape[1:])
        indices = self._get_reduce_indices()
        return reduceat(values, indices)

    def nanmin(self, values):
        return self.reduce(np.fmin.reduceat, values)

    def sum(self, values):
        return self.reduce(np.add.reduceat, values)

    def first(self, values):
        return self.reduce(lambda x, y: x[y], values)
