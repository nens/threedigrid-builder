import pygeos
from .array import array_of
import numpy as np

__all__ = ["LineStrings", "PointsOnLine"]

class LineString:
    id: int
    the_geom: pygeos.Geometry  # LineString

    
class PointOnLine:
    """A point determined by its position 'along' a linestring"""
    id: int  # just a unique number, mostly unused
    content_pk: int  # externally determined id (e.g. cs location id)
    s1d: float  # the position along the linestring
    s1d_cum: float
    linestring_idx: int


@array_of(LineString)
class LineStrings:
    @property
    def length(self):
        return pygeos.length(self.the_geom)
    
    @property
    def cum_length(self):
        result = np.zeros(len(self) + 1, dtype=float)
        np.cumsum(self.length + 1e-6, out=result[1:])
        return result

    def locate_points(self, geoms, ids, linestring_ids) -> "PointsOnLine":
        return PointsOnLine.from_geometries(self, geoms, ids, linestring_ids)
    
    def segmentize(self, desired_segment_size):
        """Divide linestrings into segments of equal length.

        Args:
            desired_segment_size (ndarray of float): the desired size of the segments; the
            actual size will depend on the linestring length and is computed by rounding
            ``line length / size`` to the nearest integer.
            Inf inputs will lead to no segmentation for that line (n_segments=1)

        Returns:
            PointsOnLine (excluding linestrings start and end)
        """
        # compute number of nodes to add per channel
        length = self.length
        n_segments = np.maximum(np.round(length / desired_segment_size).astype(int), 1)
        segment_size = length / n_segments
        n_nodes = n_segments - 1

        # get the distance to the start of each channel
        i = counts_to_row_index(n_nodes)  # e.g. [0, 0, 0, 1, 1, 3]
        j = counts_to_column_index(n_nodes)  # e.g. [0, 1, 2, 0, 1, 0]
        dist_to_start = (j + 1) * segment_size[i]

        return PointsOnLine.from_s1d(self, dist_to_start, i)


@array_of(PointOnLine)
class PointsOnLine:
    @classmethod
    def from_geometries(cls, linestrings: LineStrings, points, linestring_idx, **kwargs):
        s1d = pygeos.line_locate_point(linestrings.the_geom[linestring_idx], points)
        return cls.from_s1d(linestrings, s1d, linestring_idx, **kwargs)
    
    @classmethod
    def from_s1d(cls, linestrings: LineStrings, s1d, linestring_idx, **kwargs):
        if np.any(~np.isfinite(s1d)):
            raise ValueError("NaN values encountered in s1d")
        result = cls(
            id=np.arange(len(linestring_idx)),
            s1d=s1d,
            s1d_cum=s1d + linestrings.cum_length[linestring_idx],
            linestring_idx=linestring_idx,
            **kwargs,
        )
        result.reorder_by("s1d_cum")
        return result
    
    def as_geometries(self, line_geoms):
        return pygeos.line_interpolate_point(
            line_geoms[self.linestring_idx],
            self.s1d
        )
    
    def neighbors(self, other: "PointOnLine"):
        """Return the (indices of) the points that are before and after self.
                
        'other' must contain at least 1 point per linestring.

        If a point does not have a neigbour before it, the two returned indices
        will be equal (both pointing to the neighbour after it). Vice versa, if it
        does not have a neighboar after it, the two returned indices will be that
        of the neighboar before.
        """
        idx_2 = np.searchsorted(other.s1d_cum, self.s1d_cum)
        idx_1 = idx_2 - 1
        out_of_bounds_1 = idx_1 < 0
        out_of_bounds_2 = idx_2 >= len(other)
        idx_1[out_of_bounds_1] = 0
        idx_2[out_of_bounds_2] = len(other) - 1

        # Fix situations where idx_1 is incorrect
        extrap_mask = (other.linestring_idx[idx_1] != self.linestring_idx) | out_of_bounds_1
        idx_1[extrap_mask] = idx_2[extrap_mask]
        idx_2[extrap_mask] = np.clip(idx_2[extrap_mask] + 1, None, len(other) - 1)
        equalize = extrap_mask & (other.linestring_idx[idx_2] != self.linestring_idx)
        idx_2[equalize] -= 1

        # Fix situations where idx_2 is incorrect
        extrap_mask = (other.linestring_idx[idx_2] != self.linestring_idx) | out_of_bounds_2
        idx_2[extrap_mask] = idx_1[extrap_mask]
        idx_1[extrap_mask] = np.clip(idx_2[extrap_mask] - 1, 0, None)
        equalize = extrap_mask & (other.linestring_idx[idx_1] != self.linestring_idx)
        idx_1[equalize] += 1

        return idx_1, idx_2


def counts_to_ranges(counts):
    """Convert an array of list-of-lists counts to ranges.

    ``counts`` define lengths of lists that are concatenated into a 1D array.
    The output ranges index into that array. For example, ``arr[start[i]:stop[i]]`` will
    give all values belonging to list `i`.

    Args:
        counts (ndarray of int): list lengths

    Returns:
        tuple of start (ndarray of int), stop (ndarray of int)

    Example:
    >>> counts_to_ranges([3, 2, 0, 1, 0])
    (array([0, 3, 5, 5, 6]), array([3, 5, 5, 6, 6]))
    """
    if len(counts) == 0:
        return np.empty((0,), dtype=int), np.empty((0,), dtype=int)
    stop = np.cumsum(counts)
    start = np.roll(stop, 1)
    start[0] = 0
    return start, stop


def counts_to_row_index(counts):
    """Convert an array of list-of-lists counts into row indices into a 2D array.

    ``counts`` define lengths of lists that are concatenated into a 1D array. The output
    of this function assigns a row index to every element in the 1D array. The row index
    is counting over the outer list.

    Args:
        counts (ndarray of int): list lengths

    Returns:
        row index (ndarray of int)

    Example:
    >>> counts_to_row_index([3, 2, 0, 1, 0])
    array([0, 0, 0, 1, 1, 3])
    """
    (n,) = counts.shape
    return np.repeat(np.arange(n), counts)


def counts_to_column_index(counts):
    """Convert an array of list-of-lists counts into row indices into a 2D array.

    ``counts`` define lengths of lists that are concatenated into a 1D array. The output
    of this function assigns a column index to every element in the 1D array. The column
    index is counting over elements in a single list.

    Args:
        counts (ndarray of int): list lengths

    Returns:
        column index (ndarray of int)

    Example:
    >>> counts_to_indices([3, 2, 0, 1, 0])
    array([0, 1, 2, 0, 1, 0])
    """
    if len(counts) == 0:
        return np.empty((0,), dtype=int)
    start, stop = counts_to_ranges(counts)
    return np.arange(stop[-1]) - np.repeat(start, counts)
