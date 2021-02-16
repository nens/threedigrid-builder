import numpy as np
import pygeos


__all__ = ["Channels"]


class Channels:
    def __init__(
        self,
        the_geom,
        dist_calc_points,
        id,
        code,
        connection_node_start_id,
        connection_node_end_id,
        calculation_type,
    ):
        self.the_geom = the_geom
        self.dist_calc_points = dist_calc_points
        self.id = id
        self.code = code
        self.connection_node_start_id = connection_node_start_id
        self.connection_node_end_id = connection_node_end_id
        self.calculation_type = calculation_type

    def interpolate_nodes(self, global_dist_calc_points):
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
        dists = self.dist_calc_points.copy()  # copy because of inplace edits

        # insert default dist_calc_points where necessary
        dists[~np.isfinite(dists)] = global_dist_calc_points
        dists[dists <= 0] = global_dist_calc_points

        # compute number of nodes to add per channel
        length = pygeos.length(self.the_geom)
        n_segments = np.maximum(np.round(length / dists).astype(int), 1)
        segment_size = length / n_segments
        n_nodes = n_segments - 1
        idx = np.repeat(np.arange(self.the_geom.size), n_nodes)

        # some numpy juggling to get the distance to the start of each channel
        dist_to_start = np.arange(idx.size)
        dist_to_start[n_nodes[0] :] -= np.repeat(np.cumsum(n_nodes)[:-1], n_nodes[1:])
        dist_to_start = (dist_to_start + 1) * segment_size[idx]

        points = pygeos.line_interpolate_point(
            self.the_geom[idx],
            dist_to_start,  # note: this only copies geometry pointers
        )

        # TODO Return a to-be-implemented "Nodes" instance
        return {
            "geometry": points,
            "calculation_type": self.calculation_type[idx],
            "channel_id": self.id[idx],
            "channel_code": self.code[idx],
            "connection_node_start_id": self.connection_node_start_id[idx],
            "connection_node_end_id": self.connection_node_end_id[idx],
        }
