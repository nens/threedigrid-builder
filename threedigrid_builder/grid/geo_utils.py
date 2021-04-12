import numpy as np
import pygeos

__all__ = ["segmentize"]


def segmentize(linestrings, line_segment_size):
    """Divide linestrings into segments of equal length.

    Args:
        linestrings (ndarray of pygeos.Geometry): linestrings to segmentize
        line_segment_size (ndarray of int): the desired size of the segments. the actual size
           will depend on the total linestring length. ``length / size`` is rounded.

    Returns:
        linestrings (ndarray of pygeos.Geometry): the segments (linestrings)
        points (ndarray of pygeos.Geometry): the points where segments connect.
          this excludes the start and end of the input linestrings.
        linestring_index (ndarray of int): indices mapping segments to input linestrings
        point_index (ndarray of int): indices mapping points to input linestrings
        line_segment_size (ndarray of float): the actual length of the segments
          per input linestring
    """
    n_lines = linestrings.shape[0]
    coords, inds = pygeos.get_coordinates(linestrings, return_index=True)

    # compute the length of each vertex-vertex distance
    coord_ds = np.sqrt(np.sum((coords - np.roll(coords, 1, axis=0)) ** 2, axis=1))
    coord_s = np.cumsum(coord_ds)
    line_last_coord_idx = np.cumsum(pygeos.get_num_coordinates(linestrings)) - 1
    line_first_coord_idx = np.roll(line_last_coord_idx + 1, 1)
    line_first_coord_idx[0] = 0
    line_s = coord_s[line_last_coord_idx] - coord_s[line_first_coord_idx]

    # compute number of segments (and number of nodes) per line
    line_n_segments = np.maximum(np.round(line_s / line_segment_size).astype(int), 1)
    line_segment_size = line_s / line_n_segments
    line_n_nodes = line_n_segments - 1
    node_line_idx = np.repeat(np.arange(linestrings.size), line_n_nodes)

    # compute the "s" of each node (this matches the coord_s array)
    # first create an integer array the starts counting at 0 for every line
    _idx = np.arange(node_line_idx.size)
    _idx[line_n_nodes[0] :] -= np.repeat(np.cumsum(line_n_nodes)[:-1], line_n_nodes[1:])
    # from that, we can compute the s per line
    node_s = (_idx + 1) * line_segment_size[node_line_idx]
    node_s += np.repeat(coord_s[line_first_coord_idx], line_n_nodes)

    # interpolate the points along the linestrings
    node_insert_at = np.searchsorted(coord_s, node_s)
    node_neq_coord = ~np.isclose(coord_s[node_insert_at], node_s)
    node_add_coords = np.empty((np.count_nonzero(node_neq_coord), 2))
    node_add_coords[:, 0] = np.interp(node_s[node_neq_coord], coord_s, coords[:, 0])
    node_add_coords[:, 1] = np.interp(node_s[node_neq_coord], coord_s, coords[:, 1])

    node_coords = np.empty((node_s.size, 2))
    node_coords[~node_neq_coord] = coords[node_insert_at[~node_neq_coord]]
    node_coords[node_neq_coord, 0] = np.interp(node_s[node_neq_coord], coord_s, coords[:, 0])
    node_coords[node_neq_coord, 1] = np.interp(node_s[node_neq_coord], coord_s, coords[:, 1])

    # insert the nodes into the coordinates to generate line segments
    coords_all = np.insert(coords, node_insert_at[node_neq_coord], node_add_coords, axis=0)
    n_added_nodes = line_n_nodes - np.bincount(node_line_idx[~node_neq_coord], minlength=n_lines)
    node_neq_coord = np.roll(~node_neq_coord, 1)
    node_neq_coord[0] = False
    node_line_idx = node_insert_at + np.arange(len(node_insert_at)) - np.cumsum(node_neq_coord)
    line_last_coord_idx = line_last_coord_idx + np.cumsum(n_added_nodes)
    line_first_coord_idx = np.roll(line_last_coord_idx + 1, 1)
    line_first_coord_idx[0] = 0

    line_start_idx = np.sort(np.concatenate([line_first_coord_idx, node_line_idx]))
    line_end_idx = np.sort(np.concatenate([line_last_coord_idx, node_line_idx])) + 1
    line_coord_count = line_end_idx - line_start_idx
    line_idx = np.repeat(np.arange(len(line_coord_count)), line_coord_count)
    line_coords = coords_all[np.arange(line_coord_count.sum()) - line_idx]
    line_geometries = pygeos.linestrings(line_coords, indices=line_idx)
    node_geometries = pygeos.points(coords_all[node_line_idx])

    line_idx = np.repeat(np.arange(n_lines), line_n_nodes + 1)
    return line_geometries, node_geometries, line_idx, node_line_idx, line_segment_size
