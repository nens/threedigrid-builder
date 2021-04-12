import numpy as np
import pygeos


__all__ = ["segmentize"]


COORD_EQUAL_ATOL = 1e-8  # the distance below which coordinates are considered equal


def segmentize(linestrings, line_segment_size, return_segments=True):
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
    n_lines = linestrings.shape[0]
    coords, inds = pygeos.get_coordinates(linestrings, return_index=True)
    line_n_coords = pygeos.get_num_coordinates(linestrings)

    ## compute the length of each vertex-vertex distance
    coord_ds = np.sqrt(np.sum((coords - np.roll(coords, 1, axis=0)) ** 2, axis=1))
    coord_s = np.cumsum(coord_ds)
    line_last_coord_idx = np.cumsum(line_n_coords) - 1
    line_first_coord_idx = np.roll(line_last_coord_idx + 1, 1)
    line_first_coord_idx[0] = 0

    ## compute number of segments (and number of nodes) per line
    line_s = coord_s[line_last_coord_idx] - coord_s[line_first_coord_idx]
    # cast to single precision before rounding, the imprecision in np.sqrt may give
    # surprising results here.
    line_n_segments = np.maximum(
        np.round(line_s.astype(np.float32) / line_segment_size).astype(int), 1
    )
    line_segment_size = line_s / line_n_segments
    line_n_nodes = line_n_segments - 1
    node_line_idx = np.repeat(np.arange(linestrings.size), line_n_nodes)

    ## compute the "s" of each node (this matches the coord_s array)
    # first create an integer array the starts counting at 0 for every line
    _idx = np.arange(node_line_idx.size)
    _idx[line_n_nodes[0] :] -= np.repeat(np.cumsum(line_n_nodes)[:-1], line_n_nodes[1:])
    # from that, we can compute the s per line
    node_s = (_idx + 1) * line_segment_size[node_line_idx]
    node_s += np.repeat(coord_s[line_first_coord_idx], line_n_nodes)

    ## interpolate the points along the linestrings
    node_insert_at = np.searchsorted(coord_s, node_s)
    # find nodes that would be located exactly at an existing coordinate
    node_eq_coord = np.abs(coord_s[node_insert_at] - node_s) < COORD_EQUAL_ATOL
    node_neq_coord = ~node_eq_coord
    node_coords = np.empty((len(node_s), 2))
    node_coords[node_eq_coord] = coords[node_insert_at[node_eq_coord]]
    node_coords[node_neq_coord, 0] = np.interp(
        node_s[node_neq_coord], coord_s, coords[:, 0]
    )
    node_coords[node_neq_coord, 1] = np.interp(
        node_s[node_neq_coord], coord_s, coords[:, 1]
    )
    nodes = pygeos.points(node_coords)

    if not return_segments:
        return nodes, node_line_idx, line_segment_size

    ## insert the additional node_coords into the line coords
    segment_coords = np.insert(
        coords, node_insert_at[node_neq_coord], node_coords[node_neq_coord], axis=0
    )
    # in the line segments, node coords need to be duplicated
    _temp = np.roll(node_eq_coord, 1)
    try:
        _temp[0] = False
    except IndexError:
        pass  # edge case: no nodes
    _dup = node_insert_at + np.arange(len(node_insert_at)) - np.cumsum(_temp)
    segment_coords = np.insert(segment_coords, _dup, segment_coords[_dup], axis=0)
    # keep track where these nodes are (they are duplicated, so they are also at + 1)
    node_segment_coord_idx = _dup + np.arange(len(_dup))

    ## now find how many coordinates ended up in each segment
    # the number of added coords per input line (corrected for node-coord duplicates)
    line_n_extra_coords = line_n_nodes - np.bincount(
        node_line_idx[node_eq_coord], minlength=n_lines
    )
    # the number of coordinates in the output segments (with duplicates at each node)
    line_n_segment_coords = line_n_coords + line_n_extra_coords + line_n_nodes
    # also keep track where the line endings are
    line_last_segment_coord_idx = np.cumsum(line_n_segment_coords) - 1
    # and the line beginnings
    line_first_segment_coord_idx = np.roll(line_last_segment_coord_idx + 1, 1)
    line_first_segment_coord_idx[0] = 0
    # merge the node idx and the line first/last coordinates to get segment first/last
    segment_first_coord_idx = np.sort(
        np.concatenate([line_first_segment_coord_idx, node_segment_coord_idx + 1])
    )
    segment_last_coord_idx = np.sort(
        np.concatenate([line_last_segment_coord_idx, node_segment_coord_idx])
    )
    # and finally the number of coordinates per segment
    segment_n_coords = segment_last_coord_idx - segment_first_coord_idx + 1

    ## construct the segments
    segments = pygeos.linestrings(
        segment_coords,
        indices=np.repeat(np.arange(segment_n_coords.shape[0]), segment_n_coords),
    )

    ## find which line each segment belongs to
    segment_line_idx = np.repeat(np.arange(len(line_n_segments)), line_n_segments)

    return nodes, node_line_idx, line_segment_size, segments, segment_line_idx
