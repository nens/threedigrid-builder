import pygeos
from .array import Array
import numpy as np

__all__ = ["LineStrings", "PointsOnLine", "LinesOnLine"]


COORD_EQUAL_ATOL = 1e-8  # the distance below which coordinates are considered equal

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

class LineOnLine:
    """A line determined by its position 'along' a linestring"""
    id: int  # just a unique number, mostly unused
    s1d_start: float  # the position along the linestring
    s1d_end: float  # the position along the linestring
    linestring_idx: int

class LineStrings(Array[LineString]):
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
    
    def interpolate_points(self, desired_segment_size):
        """Compute points that divide linestrings into segments of equal length.

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
    
    def segmentize(self, points: "PointsOnLine") -> "LineOnLine":
        """Return lines that result from splitting self at 'points'
        """
        segment_counts = np.bincount(points.linestring_idx, minlength=len(self)) + 1

        # cut the channel geometries into segment geometries
        start_s, end_s, segment_idx = segment_start_end(
            self.the_geom, segment_counts, points.s1d
        )
        return LinesOnLine(id=range(len(segment_idx)), s1d_start=start_s, s1d_end=end_s, linestring_idx=segment_idx)


class PointsOnLine(Array[PointOnLine]):
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

class LinesOnLine(Array[LineOnLine]):
    @property
    def s1d(self):
        return (self.s1d_start + self.s1d_end) / 2
    
    @property
    def ds1d(self):
        return self.s1d_end - self.s1d_start

    @property
    def ds1d_half(self):
        return self.ds1d / 2

    def as_geometries(self, line_geoms):
        return line_substring(line_geoms, self.s1d_start, self.s1d_end, self.linestring_idx)


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



def line_substring(linestrings, start, end, index=None):
    """Divide linestrings into segments given [start, end] measured along the line.

    Note that there is no bound check done on start and end. Expect bogus results or
    errors if start is negative or if end is larger than the linestring length.

    Args:
        linestrings (ndarray of pygeos.Geometry): Linestrings to segmentize
        start (ndarray of float): The start of the segment, measured along the line.
        end (ndarray of float): The end of the segment, measured along the line.
        index (ndarray of int, optional): An optional linestring index per start/end.
           If not given, then all input arrays should be equal sized.

    Returns:
        segments: the segments (sublinestrings of the input linestrings)
    """
    (n_lines,) = linestrings.shape
    if index is None:
        n_segments = n_lines
        index = np.arange(n_lines)
    else:
        (n_segments,) = index.shape

    if n_segments == 0:
        return np.empty((0,), dtype=object)

    coords, coord_line_idx = pygeos.get_coordinates(linestrings, return_index=True)
    line_n_coords = pygeos.get_num_coordinates(linestrings)

    ## compute the length of each vertex-vertex distance
    coord_ds = np.sqrt(np.sum((coords - np.roll(coords, 1, axis=0)) ** 2, axis=1))
    coord_s = np.cumsum(coord_ds)
    # compute corresponding length of each given start and end
    line_first_coord_idx, line_end_coord_idx = counts_to_ranges(line_n_coords)
    start_s = start + coord_s[line_first_coord_idx][index]
    end_s = end + coord_s[line_first_coord_idx][index]

    ## interpolate the start & end points along the linestrings
    start_insert_at = np.searchsorted(coord_s, start_s - COORD_EQUAL_ATOL)
    end_insert_at = np.searchsorted(coord_s, end_s - COORD_EQUAL_ATOL)
    # find points that would be located exactly at an existing coordinate
    start_eq_coord = np.abs(coord_s[start_insert_at] - start_s) < COORD_EQUAL_ATOL
    end_eq_coord = np.abs(coord_s[end_insert_at] - end_s) < COORD_EQUAL_ATOL
    # insert the existing coords and interpolate the others
    start_coords = np.empty((n_segments, 2))
    end_coords = np.empty((n_segments, 2))
    start_coords[start_eq_coord] = coords[start_insert_at[start_eq_coord]]
    end_coords[end_eq_coord] = coords[end_insert_at[end_eq_coord]]
    for i in range(2):
        start_coords[~start_eq_coord, i] = np.interp(
            start_s[~start_eq_coord], coord_s, coords[:, i]
        )
        end_coords[~end_eq_coord, i] = np.interp(
            end_s[~end_eq_coord], coord_s, coords[:, i]
        )

    ## fill an array of coordinates for the segmetns
    # the number of coordinates per segment are:
    # + 2 (line start and end)
    # + coordinates in between (end_insert_at - start_insert_at)
    # + correction where start_eq_coord (then there is 1 less "in between")
    # + no correction where end_eq_coord
    segment_n_coords = 2 + end_insert_at - start_insert_at - start_eq_coord
    segment_first_coord_idx, segment_end_coord_idx = counts_to_ranges(segment_n_coords)

    # fill in the start & end coords
    segment_coords = np.empty((segment_n_coords.sum(), 2))
    segment_coords[segment_first_coord_idx] = start_coords
    segment_coords[segment_end_coord_idx - 1] = end_coords

    # fill the remaining segment coordinates with original coordinates
    i = counts_to_row_index(segment_n_coords - 2)  # e.g. [0, 0, 0, 1, 1, 3]
    j = counts_to_column_index(segment_n_coords - 2)  # e.g. [0, 1, 2, 0, 1, 0]
    segment_coords_to_fill = segment_first_coord_idx[i] + j + 1
    coords_to_fill_with = start_insert_at[i] + start_eq_coord[i] + j
    segment_coords[segment_coords_to_fill] = coords[coords_to_fill_with]

    # construct the segments
    segments = pygeos.linestrings(
        segment_coords,
        indices=counts_to_row_index(segment_n_coords),
    )

    return segments


def segment_start_end(linestrings, segment_counts, dist_to_start):
    """Utility function to use the output of segmentize as input of line_substrings.

    Args:
        linestrings (ndarray of pygeos.Geometry): linestrings to segmentize
        segment_counts (ndarray of int): the number of segments per linestring
        dist_to_start (ndarray of float): the location of added nodes, measured along
          the linestring. The length of this array is such that:
          ``len(linestrings) + len(dist_to_start) = segment_counts.sum()``

    Returns:
        tuple of:
        - segment_start (ndarray of float): The segment start, measured along the line.
        - segment_end (ndarray of float): The segment end, measured along the line.
        - segment_index (ndarray of int): Indices mapping segments to input linestrings.
    """
    start_idx, last_idx = counts_to_ranges(segment_counts)
    last_idx -= 1

    start_s = np.full(segment_counts.sum(), np.nan)
    end_s = np.full(segment_counts.sum(), np.nan)
    # every first start_s is 0.0, the rest are added nodes with an s1d
    start_s[start_idx] = 0.0
    start_s[np.isnan(start_s)] = dist_to_start
    # every last end_s equals to the length, the rest are added nodes with an s1d
    end_s[last_idx] = pygeos.length(linestrings)
    end_s[np.isnan(end_s)] = dist_to_start
    return start_s, end_s, counts_to_row_index(segment_counts)
