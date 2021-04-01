import numpy as np
import pygeos


def segmentize(linestrings, segment_size):
    """Divide linestrings into segments of equal length.

    Args:
        linestrings (ndarray of pygeos.Geometry): linestrings to segmentize
        segment_size (ndarray of int): the desired size of the segments. the actual size
           will depend on the total linestring length. ``length / size`` is rounded.

    Returns:
        linestrings (ndarray of pygeos.Geometry): the segments (linestrings)
        points (ndarray of pygeos.Geometry): the points that divide the input linestrings
        linestring_index (ndarray of int): indices mapping segments to input linestrings
        point_index (ndarray of int): indices mapping points to input linestrings
        segment_size (ndarray of float): the actual length of the segments
          per input linestring
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
    node_neq_coord = ~np.isclose(s[insert_at], s_nodes)
    coords_nodes = np.empty((np.count_nonzero(node_neq_coord), 2))
    coords_nodes[:, 0] = np.interp(s_nodes[node_neq_coord], s, coords[:, 0])
    coords_nodes[:, 1] = np.interp(s_nodes[node_neq_coord], s, coords[:, 1])

    # insert the nodes into the coordinates to generate line segments
    coords_all = np.insert(coords, insert_at[node_neq_coord], coords_nodes, axis=0)
    n_added_nodes = n_nodes - np.bincount(idx[~node_neq_coord], minlength=n_lines)
    node_neq_coord = np.roll(~node_neq_coord, 1)
    node_neq_coord[0] = False
    node_idx = insert_at + np.arange(len(insert_at)) - np.cumsum(node_neq_coord)
    end_idx = end_idx + np.cumsum(n_added_nodes)
    start_idx = np.roll(end_idx + 1, 1)
    start_idx[0] = 0

    line_start_idx = np.sort(np.concatenate([start_idx, node_idx]))
    line_end_idx = np.sort(np.concatenate([end_idx, node_idx])) + 1
    line_coord_count = line_end_idx - line_start_idx
    line_idx = np.repeat(np.arange(len(line_coord_count)), line_coord_count)
    line_coords = coords_all[np.arange(line_coord_count.sum()) - line_idx]
    line_geometries = pygeos.linestrings(line_coords, indices=line_idx)
    node_geometries = pygeos.points(coords_all[node_idx])

    line_idx = np.repeat(np.arange(n_lines), n_nodes + 1)
    return line_geometries, node_geometries, line_idx, idx, segment_size
