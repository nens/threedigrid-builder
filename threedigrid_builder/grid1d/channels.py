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

    def __repr__(self):
        return "<Channels object (len:{})>".format(len(self.the_geom))

    def interpolate_nodes(self, global_dist_calc_points):
        """Compute interpolated channel nodes

        Args:
          global_dist_calc_points (float): Default node interdistance.

        Returns:
          dict of nodes with the following properties (all 1D arrays):
          - geometry
          - calculation_type
          - channel_id
          - channel_code
          - connection_node_start_id
          - connection_node_end_id
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

    def nodes_and_lines(
        self, global_dist_calc_points, channel_node_offset, connection_node_offset
    ):
        """Compute all channel nodes and lines

        Args:
          connection_nodes (ConnectionNodes): all available connection nodes
          global_dist_calc_points (float): Default node interdistance.
          channel_node_offset (int): the index of the first channel node in the
            target node array (for the lines)
          connection_node_offset (int): the index of the first connection node in the
            target node array (for the lines)

        Returns:
          dict of nodes with the following properties (all 1D arrays):
          - geometry
          - calculation_type
          - channel_id
          - channel_code
          - connection_node_start_id
          - connection_node_end_id
          dict of lines with the following properties:
          - nodes (N, 2) array of node ids
        """
        nodes = self.interpolate_nodes(global_dist_calc_points)
        n_nodes = len(nodes["geometry"])
        lines = np.vstack(
            [np.arange(n_nodes), np.arange(1, n_nodes + 1)]
        ) + channel_node_offset

        # ending nodes should not connect to the next channel, but to the
        # connection_node_end_id of that channel
        is_channel_end = np.append(
            nodes["channel_id"][1:] != nodes["channel_id"][:-1], [True]
        )
        lines[1][is_channel_end] = nodes["connection_node_end_id"][is_channel_end] + connection_node_offset


        # Todo add starting lines (for channels without and with interpolated node)
