import numpy as np
import pygeos


__all__ = ["segmentize", "line_substring"]


COORD_EQUAL_ATOL = 1e-8  # the distance below which coordinates are considered equal


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
    start, stop = counts_to_ranges(counts)
    return np.arange(stop[-1]) - np.repeat(start, counts)


def segmentize(linestrings, desired_segment_size):
    """Return points that divide linestrings into segments of equal length.

    Args:
        linestrings (ndarray of pygeos.Geometry): linestrings to segmentize
        desired_segment_size (ndarray of float): the desired size of the segments; the
           actual size will depend on the linestring length and is computed by rounding
           ``line length / size`` to the nearest integer.

    Returns:
        nodes: the points where segments connect.
          this excludes the start and end of the input linestrings.
        line_idx: indices mapping nodes to input linestrings
        segment_size: the actual length of the segments per input linestring
    """
    assert linestrings.ndim == desired_segment_size.ndim == 1
    assert linestrings.shape == desired_segment_size.shape

    # compute number of nodes to add per channel
    length = pygeos.length(linestrings)
    n_segments = np.maximum(np.round(length / desired_segment_size).astype(int), 1)
    segment_size = length / n_segments
    n_nodes = n_segments - 1

    # get the distance to the start of each channel
    i = counts_to_row_index(n_nodes)  # e.g. [0, 0, 0, 1, 1, 3]
    j = counts_to_column_index(n_nodes)  # e.g. [0, 1, 2, 0, 1, 0]
    dist_to_start = (j + 1) * segment_size[i]

    nodes = pygeos.line_interpolate_point(
        linestrings[i],
        dist_to_start,  # note: this only copies geometry pointers
    )
    return nodes, i


def line_substring(linestrings, start, end, index=None):
    """Divide linestrings into segments given [start, end] measured along the line.

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
        assert start.shape == end.shape == (n_lines,)
        index = np.arange(n_lines)
    else:
        (n_segments,) = index.shape
        assert start.shape == end.shape == (n_segments,)

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
    start_insert_at = np.searchsorted(coord_s, start_s - 1e-15)
    end_insert_at = np.searchsorted(coord_s, end_s - 1e-15)
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


def segment_start_end(linestrings, segment_counts):
    """Utility function to use the output of segmentize as input of line_substrings.

    Args:
        linestrings (ndarray of pygeos.Geometry): linestrings to segmentize
        segment_counts (ndarray of int): the number of segments per linestring

    Returns:
        tuple of:
        - segment_start (ndarray of float): The segment start, measured along the line.
        - segment_end (ndarray of float): The segment end, measured along the line.
        - segment_index (ndarray of int): Indices mapping segments to input linestrings.
    """
    assert linestrings.ndim == segment_counts.ndim == 1
    assert linestrings.shape == segment_counts.shape
    i = counts_to_row_index(segment_counts)  # e.g. [0, 0, 0, 1, 1, 3]
    j = counts_to_column_index(segment_counts)  # e.g. [0, 1, 2, 0, 1, 0]

    lengths = pygeos.length(linestrings)
    segment_size = lengths / segment_counts
    mapped_segment_size = (segment_size)[i]
    start = j * mapped_segment_size
    end = start + mapped_segment_size
    return start, end, i
