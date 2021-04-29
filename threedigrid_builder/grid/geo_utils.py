import numpy as np
import pygeos


__all__ = ["segmentize", "line_substring"]


COORD_EQUAL_ATOL = 1e-8  # the distance below which coordinates are considered equal


def segmentize(linestrings, line_segment_size):
    """Return points that divide linestrings into segments of equal length.

    Args:
        linestrings (ndarray of pygeos.Geometry): linestrings to segmentize
        line_segment_size (ndarray of int): the desired size of the segments. the actual size
           will depend on the total linestring length. ``length / size`` is rounded.

    Returns:
        nodes: the points where segments connect.
          this excludes the start and end of the input linestrings.
        line_idx: indices mapping nodes to input linestrings
        segment_size: the actual length of the segments per input linestring
    """
    n_linestrings = linestrings.size

    # compute number of nodes to add per channel
    length = pygeos.length(linestrings)
    n_segments = np.maximum(np.round(length / line_segment_size).astype(int), 1)
    segment_size = length / n_segments
    n_nodes = n_segments - 1
    line_idx = np.repeat(np.arange(n_linestrings), n_nodes)

    # some numpy juggling to get the distance to the start of each channel
    dist_to_start = np.arange(line_idx.size)
    dist_to_start[n_nodes[0] :] -= np.repeat(np.cumsum(n_nodes)[:-1], n_nodes[1:])
    dist_to_start = (dist_to_start + 1) * segment_size[line_idx]

    nodes = pygeos.line_interpolate_point(
        linestrings[line_idx],
        dist_to_start,  # note: this only copies geometry pointers
    )
    return nodes, line_idx, segment_size


def line_substring(linestrings, start, end, index=None):
    """Divide linestrings into segments.

    Args:
        linestrings (ndarray of pygeos.Geometry): linestrings to segmentize
        start (ndarray of float): the start of the segment, measured along the line
        end (ndarray of float): the end of the segment, measured along the line
        index (ndarray of int, optional): an optional mapping giving the linestring
           index per start/end. If not given, linestrings and start/end should be equal
           sized.

    Returns:
        segments: the segments
    """
    n_lines, = linestrings.shape
    if index is None:
        n_segments = n_lines
        assert start.shape == end.shape == (n_lines, )
        index = np.arange(n_lines)
    else:
        n_segments, = index.shape
        assert start.shape == end.shape == (n_segments, )

    coords, coord_line_idx = pygeos.get_coordinates(linestrings, return_index=True)
    line_n_coords = pygeos.get_num_coordinates(linestrings)

    ## compute the length of each vertex-vertex distance
    coord_ds = np.sqrt(np.sum((coords - np.roll(coords, 1, axis=0)) ** 2, axis=1))
    coord_s = np.cumsum(coord_ds)
    line_last_coord_idx = np.cumsum(line_n_coords) - 1
    line_first_coord_idx = np.roll(line_last_coord_idx + 1, 1)
    line_first_coord_idx[0] = 0
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

    # the coordinates per segment are:
    # + 2 (line start and end)
    # + coordinates in between (end_insert_at - start_insert_at)
    # + correction where start_eq_coord (then there is 1 less "in between")
    # + no correction where end_eq_coord
    segment_n_coords = 2 + end_insert_at - start_insert_at - start_eq_coord
    segment_coords = np.empty((segment_n_coords.sum(), 2))
    segment_last_coord_idx = np.cumsum(segment_n_coords) - 1
    segment_first_coord_idx = np.roll(segment_last_coord_idx + 1, 1)
    segment_first_coord_idx[0] = 0

    segment_coords[segment_first_coord_idx] = start_coords
    segment_coords[segment_last_coord_idx] = end_coords

    # generate indices with i = segment idx and j = coord within segment idx
    # e.g. i=[0, 0, 0, 0, 1, 2, 2]; j=[0, 1, 2, 3, 0, 0, 1] for counts = [4, 1, 2]
    _first_idx = np.roll(np.cumsum(segment_n_coords - 2), 1)
    _total, _first_idx[0] = _first_idx[0], 0
    i = np.repeat(np.arange(n_segments), segment_n_coords - 2)
    j = np.arange(_total) - np.repeat(_first_idx, segment_n_coords - 2)

    # fill the remaining segment coordinates with original coordinates
    segment_coords_to_fill = segment_first_coord_idx[i] + j + 1
    coords_to_fill_with = start_insert_at[i] + start_eq_coord[i] + j
    segment_coords[segment_coords_to_fill] = coords[coords_to_fill_with]

    ## construct the segments
    segments = pygeos.linestrings(
        segment_coords,
        indices=np.repeat(np.arange(n_segments), segment_n_coords),
    )

    return segments


def segment_size_to_start_end(linestrings, node_line_idx, line_segment_size):
    """Divide linestrings into segments of equal length.

    Optionally returns the segments that, pieced together, create the linestrings.

    Args:
        linestrings (ndarray of pygeos.Geometry): linestrings to segmentize
        line_segment_size (ndarray of int): the desired size of the segments. the actual size
           will depend on the total linestring length. ``length / size`` is rounded.
        return_segments (bool): whether to return the segments, default True

    Returns:
        nodes: the points where segments connect.
          this excludes the start and end of the input linestrings.
        node_index: indices mapping nodes to input linestrings
        segment_size: the actual length of the segments per input linestring

        only if return_segments is True:
            segments: the segments (linestrings)
            segment_index: indices mapping segments to input linestrings
    """
    nodes, node_line_idx, actual_line_segment_size = segmentize(linestrings, line_segment_size)

    n_lines = linestrings.shape[0]
    segment_counts = np.bincount(node_line_idx, minlength=n_lines) + 1
    segment_line_idx = np.repeat(np.arange(n_lines), segment_counts)

    # generate an array of indices that increase within every linestring
    # e.g. [0, 1, 2, 3, 0, 1, 0, 1, 2, 3, 4, 5] for segment_counts = [4, 2, 6]
    first_idx = np.roll(np.cumsum(segment_counts), 1)
    total, first_idx[0] = first_idx[0], 0
    line_segment_idx = np.arange(total) - np.repeat(first_idx, segment_counts)

    segment_start = line_segment_idx * actual_line_segment_size[segment_line_idx]
    segment_end = segment_start + actual_line_segment_size[segment_line_idx]
    return segment_start, segment_end, segment_line_idx
