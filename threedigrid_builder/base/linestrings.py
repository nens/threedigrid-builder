import numpy as np
import shapely

from .array import Array

__all__ = ["LineStrings", "PointsOnLine", "LinesOnLine"]


COORD_EQUAL_ATOL = 1e-8  # the distance below which coordinates are considered equal


class PointOnLine:
    """A point determined by its position 'along' a linestring"""

    id: int  # just a unique number, mostly unused
    content_pk: int  # externally determined id
    s1d: float  # the position along the linestring
    linestring_idx: int
    secondary_content_pk: int  # for double connected points


class LineOnLine:
    """A line determined by its position 'along' a linestring"""

    id: int  # just a unique number, mostly unused
    content_pk: int  # externally determined id
    s1d_start: float  # the position along the linestring
    s1d_end: float  # the position along the linestring
    linestring_idx: int


class LineStrings:
    def __init__(self, objects):
        geom_types = shapely.get_type_id(objects.the_geom)
        assert np.all((geom_types == -1) | (geom_types == 1))
        self.objects = objects

    def __len__(self):
        return len(self.objects)

    def __getattr__(self, name):
        return getattr(self.objects, name)

    @property
    def length(self):
        return shapely.length(self.the_geom)

    @property
    def length_cumulative(self):
        result = np.zeros(len(self) + 1, dtype=float)
        np.cumsum(self.length + 1e-6, out=result[1:])
        return result

    def segmentize(self, points: "PointsOnLine") -> "LinesOnLine":
        """Return lines that result from splitting self at 'points'"""
        segment_counts = np.bincount(points.linestring_idx, minlength=len(self)) + 1

        # cut the channel geometries into segment geometries
        start_s, end_s, segment_idx = segment_start_end(
            self.the_geom, segment_counts, points.s1d
        )
        return LinesOnLine(
            linestrings=self,
            id=range(len(segment_idx)),
            s1d_start=start_s,
            s1d_end=end_s,
            linestring_idx=segment_idx,
        )

    def sanitize(self, points_1, points_2) -> None:
        """Make sure that the self go from points_1 to points_2.

        The linestring is reversed if the distance from the linestring start to point_1
        and v.v. is lowered.

        Then, point_1 and/or point_2 are added to the linestring if the distance is larger
        than the tolerance of 1e-8.

        This function acts inplace on self.
        """
        has_geom = np.where(shapely.is_geometry(self.the_geom))[0]
        if len(has_geom) == 0:
            return

        # reverse the geometries where necessary
        node_start = points_1[has_geom]
        node_end = points_2[has_geom]
        line_start = shapely.get_point(self.the_geom[has_geom], 0)
        line_end = shapely.get_point(self.the_geom[has_geom], -1)

        dist_start_start = shapely.distance(node_start, line_start)
        dist_end_end = shapely.distance(node_end, line_end)
        dist_start_end = shapely.distance(node_start, line_end)
        dist_end_start = shapely.distance(node_end, line_start)
        needs_reversion = has_geom[
            (dist_start_start > dist_start_end) & (dist_end_end > dist_end_start)
        ]
        if len(needs_reversion) > 0:
            self.the_geom[needs_reversion] = shapely.reverse(
                self.the_geom[needs_reversion]
            )
            line_start[needs_reversion] = shapely.get_point(
                self.the_geom[needs_reversion], 0
            )
            line_end[needs_reversion] = shapely.get_point(
                self.the_geom[needs_reversion], -1
            )
            dist_start_start[needs_reversion] = dist_start_end[needs_reversion]
            dist_end_end[needs_reversion] = dist_end_start[needs_reversion]

        # for the remaining geometries; add point if necessary
        add_first = has_geom[dist_start_start > COORD_EQUAL_ATOL]
        if len(add_first) > 0:
            self.the_geom[add_first] = prepend_point(
                self.the_geom[add_first], points_1[add_first]
            )
        add_end = has_geom[dist_end_end > COORD_EQUAL_ATOL]
        if len(add_end) > 0:
            self.the_geom[add_end] = append_point(
                self.the_geom[add_end], points_2[add_end]
            )


class PointsOnLine(Array[PointOnLine]):
    scalars = ("linestrings",)

    @property
    def linestring_id(self):
        return self.linestrings.index_to_id(self.linestring_idx)

    @classmethod
    def empty(cls, linestrings: LineStrings):
        return cls(linestrings=linestrings, id=[])

    @classmethod
    def from_geometries(
        cls, linestrings: LineStrings, points, linestring_idx, **kwargs
    ):
        s1d = shapely.line_locate_point(linestrings.the_geom[linestring_idx], points)
        return cls.from_s1d(linestrings, s1d, linestring_idx, **kwargs)

    @classmethod
    def from_s1d(cls, linestrings: LineStrings, s1d, linestring_idx, **kwargs):
        if np.any(~np.isfinite(s1d)):
            raise ValueError("NaN values encountered in s1d")
        result = cls(
            linestrings=linestrings,
            id=np.arange(len(linestring_idx)),
            s1d=s1d,
            linestring_idx=linestring_idx,
            **kwargs,
        )
        result.reorder(np.argsort(result.s1d_cum))
        return result

    @property
    def s1d_cum(self):
        return self.s1d + self.linestrings.length_cumulative[self.linestring_idx]

    @property
    def the_geom(self):
        return shapely.line_interpolate_point(
            self.linestrings.the_geom[self.linestring_idx], self.s1d
        )

    @property
    def at_start(self):
        return self.s1d == 0.0

    @property
    def at_end(self):
        return self.s1d == self.linestrings.length[self.linestring_idx]

    def merge_with(self, other: "PointsOnLine"):
        if len(other) == 0:
            return self
        if len(self) == 0:
            first_i = 1
        else:
            first_i = self.id[-1] + 1
        other = other[:]
        other.id[:] = np.arange(first_i, first_i + len(other))
        result = self + other
        result.reorder(np.argsort(result.s1d_cum))
        return result

    def neighbours(self, other: "PointsOnLine"):
        """Compute neighbours in self of 'other'

        Returns two arrays that index into self: the neighbours 'before' and the
        neighbours 'after'.

        'other' must contain at least 1 point per linestring.

        If a point does not have a neigbour before it, the two returned indices
        will be the two points after it. Vice versa, if it does not have a neighbour after
        it, the two returned indices will be that of the neighboar before.

        If there there is only 1 point on a linestring, the neighbours will be equal.
        """
        idx_2 = np.searchsorted(self.s1d_cum, other.s1d_cum)
        idx_1 = idx_2 - 1
        out_of_bounds_1 = idx_1 < 0
        out_of_bounds_2 = idx_2 >= len(self)
        idx_1[out_of_bounds_1] = 0
        idx_2[out_of_bounds_2] = len(self) - 1

        # Fix situations where idx_1 is incorrect
        extrap_mask = (
            self.linestring_idx[idx_1] != other.linestring_idx
        ) | out_of_bounds_1
        idx_1[extrap_mask] = idx_2[extrap_mask]
        idx_2[extrap_mask] = np.clip(idx_2[extrap_mask] + 1, None, len(self) - 1)
        equalize = extrap_mask & (self.linestring_idx[idx_2] != other.linestring_idx)
        idx_2[equalize] -= 1

        # Fix situations where idx_2 is incorrect
        extrap_mask = (
            self.linestring_idx[idx_2] != other.linestring_idx
        ) | out_of_bounds_2
        idx_2[extrap_mask] = idx_1[extrap_mask]
        idx_1[extrap_mask] = np.clip(idx_2[extrap_mask] - 1, 0, None)
        equalize = extrap_mask & (self.linestring_idx[idx_1] != other.linestring_idx)
        idx_1[equalize] += 1

        return idx_1, idx_2


class LinesOnLine(Array[LineOnLine]):
    scalars = ("linestrings",)

    @property
    def s1d(self):
        return (self.s1d_start + self.s1d_end) / 2

    @property
    def ds1d(self):
        return self.s1d_end - self.s1d_start

    @property
    def the_geom(self):
        return line_substring(
            self.linestrings.the_geom, self.s1d_start, self.s1d_end, self.linestring_idx
        )

    def interpolate_points(self, desired_segment_size):
        """Compute points that divide sublinestrings into segments of equal length.

        Args:
            desired_segment_size (ndarray of float): the desired size of the segments; the
            actual size will depend on the linestring length and is computed by rounding
            ``line length / size`` to the nearest integer.
            Inf inputs will lead to no segmentation for that line (n_segments=1)

        Returns:
            PointsOnLine (excluding linestrings start and end)
        """
        # compute number of nodes to add per channel
        length = self.ds1d
        n_segments = np.maximum(np.round(length / desired_segment_size).astype(int), 1)
        segment_size = length / n_segments
        n_nodes = n_segments - 1

        # get the distance to the start of each channel
        i = counts_to_row_index(n_nodes)  # e.g. [0, 0, 0, 1, 1, 3]
        j = counts_to_column_index(n_nodes)  # e.g. [0, 1, 2, 0, 1, 0]
        dist_to_start = (j + 1) * segment_size[i] + self.s1d_start[i]

        return PointsOnLine.from_s1d(
            self.linestrings, dist_to_start, self.linestring_idx[i]
        )


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
        linestrings (ndarray of shapely.Geometry): Linestrings to segmentize
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

    coords, coord_line_idx = shapely.get_coordinates(linestrings, return_index=True)
    line_n_coords = shapely.get_num_coordinates(linestrings)

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
    segments = shapely.linestrings(
        segment_coords,
        indices=counts_to_row_index(segment_n_coords),
    )

    return segments


def segment_start_end(linestrings, segment_counts, dist_to_start):
    """Utility function to use the output of segmentize as input of line_substrings.

    Args:
        linestrings (ndarray of shapely.Geometry): linestrings to segmentize
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
    end_s[last_idx] = shapely.length(linestrings)
    end_s[np.isnan(end_s)] = dist_to_start
    return start_s, end_s, counts_to_row_index(segment_counts)


def _add_point(linestrings, points, at=0):
    """Append or prepend points to linestrings"""
    counts = shapely.get_num_coordinates(linestrings)
    coords, index = shapely.get_coordinates(linestrings, return_index=True)
    start, end = counts_to_ranges(counts)
    if at == 0:
        insert_at = start
    elif at == -1:
        insert_at = end
    else:
        raise ValueError(f"Cannot insert at {at}")
    new_coords = np.insert(coords, insert_at, shapely.get_coordinates(points), axis=0)
    new_index = np.insert(index, insert_at, np.arange(len(insert_at)))
    return shapely.linestrings(new_coords, indices=new_index)


def prepend_point(linestrings, points):
    return _add_point(linestrings, points, at=0)


def append_point(linestrings, points):
    return _add_point(linestrings, points, at=-1)
