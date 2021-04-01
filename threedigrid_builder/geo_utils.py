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
    n_lines = linestrings.shape[0]
    coords, inds = pygeos.get_coordinates(linestrings, return_index=True)
    n_points = pygeos.get_num_coordinates(linestrings)

    ds = np.sqrt(np.sum((coords - np.roll(coords, 1, axis=0)) ** 2, axis=1))
    s = np.cumsum(ds)
    end_idx = np.cumsum(n_points) - 1
    start_idx = np.roll(end_idx + 1, 1)
    start_idx[0] = 0

    # compute number of nodes to add per line
    length = s[end_idx] - s[start_idx]
    n_segments = np.maximum(np.round(length / segment_size).astype(int), 1)
    segment_size = length / n_segments
    n_nodes = n_segments - 1
    idx = np.repeat(np.arange(linestrings.size), n_nodes)

    # some numpy juggling to get the distance to the start of each
    dist_to_start = np.arange(idx.size)
    dist_to_start[n_nodes[0] :] -= np.repeat(np.cumsum(n_nodes)[:-1], n_nodes[1:])
    dist_to_start = (dist_to_start + 1) * segment_size[idx]

    # interpolate the points along the linestrings
    s_nodes = np.repeat(s[start_idx], n_nodes) + dist_to_start
    insert_at = np.searchsorted(s, s_nodes)
    mask = ~np.isclose(s[insert_at], s_nodes)
    coords_nodes = np.empty((np.count_nonzero(mask), 2))
    coords_nodes[:, 0] = np.interp(s_nodes[mask], s, coords[:, 0])
    coords_nodes[:, 1] = np.interp(s_nodes[mask], s, coords[:, 1])

    # make the coordinates of the lines
    2 * n_lines + 2 * n_nodes 

    coords_all = np.insert(coords, insert_at[mask], coords_nodes, axis=0)
    n_added_nodes = n_nodes - np.bincount(idx[~mask], minlength=n_lines)

    # find the indices of the node coordinates in the complete coord array
    mask2 = np.roll(~mask, 1)
    mask2[0] = False
    node_idx = insert_at + np.arange(len(insert_at)) - np.cumsum(mask2)
    channel_idx = np.repeat(np.arange(n_lines), n_points + n_added_nodes)
    end_idx = end_idx + np.cumsum(n_added_nodes)
    start_idx = np.roll(end_idx + 1, 1)
    start_idx[0] = 0

    # assemble the coordinates for each segment (this will duplicate nodes)
        

    points = pygeos.line_interpolate_point(linestrings[idx], dist_to_start)
    return points, segment_size, idx
