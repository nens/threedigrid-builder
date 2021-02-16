from threedigrid_builder.interface import get_channels
from threedigrid_builder.interface import get_global_settings

import numpy as np
import pygeos


__all__ = ["node_channels"]


def node_channels(path):
    """Compute interpolated channel nodes

    Args:
      path (str): Path to an SQLite

    Returns:
      dict with the following keys:
      - geometry (array of points)
      - dist_calc_points (array of )
      - code (ndarray of python str objects)
      - display_name (ndarray of python str objects)
      - calculation_type (ndarray of uint8)
    """
    # load data
    data = get_channels(path)
    geoms = data["the_geom"]
    dists = data["dist_calc_points"]

    # insert default dist_calc_points where necessary
    global_dists = get_global_settings(path)["dist_calc_points"]
    dists[~np.isfinite(dists)] = global_dists
    dists[dists <= 0] = global_dists

    # compute number of nodes to add per channel
    length = pygeos.length(geoms)
    n_segments = np.maximum(np.round(length / data["dist_calc_points"]).astype(int), 1)
    segment_size = length / n_segments
    n_nodes = n_segments - 1
    idx = np.repeat(np.arange(geoms.size), n_nodes)

    # some numpy juggling to get the distance to the start of each channel
    dist_to_start = np.arange(idx.size)
    dist_to_start[n_nodes[0] :] -= np.repeat(np.cumsum(n_nodes)[:-1], n_nodes[1:])
    dist_to_start = (dist_to_start + 1) * segment_size[idx]

    points = pygeos.line_interpolate_point(
        geoms[idx], dist_to_start  # note: this only copies geometry pointers
    )

    return {
        "geometry": points,
        "calculation_type": data["calculation_type"][idx],
        "channel_id": data["id"][idx],
        "channel_code": data["code"][idx],
        "connection_node_start_id": data["connection_node_start_id"][idx],
        "connection_node_end_id": data["connection_node_end_id"][idx],
    }
