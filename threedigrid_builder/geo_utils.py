import numpy as np
import pygeos


def segmentize(linestrings, segment_size):
    """Divide linestrings into segments of equal length.

    Args:
        linestrings (ndarray of pygeos.Geometry): linestrings to segmentize
        segment_size (ndarray of int): the desired size of the segments. the actual size
           will depend on the total linestring length. ``length / size`` is rounded.

    Returns:
        points (ndarray of pygeos.Geometry): array size equals the sum of n_segments
        segment_size (ndarray of float): the actual length of the segments
        index (ndarray of int): indices mapping points to linestrings
    """
    # compute number of nodes to add per line
    length = pygeos.length(linestrings)
    n_segments = np.maximum(np.round(length / segment_size).astype(int), 1)
    segment_size = length / n_segments
    n_nodes = n_segments - 1
    idx = np.repeat(np.arange(linestrings.size), n_nodes)

    # some numpy juggling to get the distance to the start of each
    dist_to_start = np.arange(idx.size)
    dist_to_start[n_nodes[0] :] -= np.repeat(np.cumsum(n_nodes)[:-1], n_nodes[1:])
    dist_to_start = (dist_to_start + 1) * segment_size[idx]

    # interpolate the points along the linestrings
    points = pygeos.line_interpolate_point(linestrings[idx], dist_to_start)
    return points, segment_size, idx
