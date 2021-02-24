from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import Grid
from threedigrid_builder.grid import Lines
from threedigrid_builder.grid import Nodes

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
            self (Channels): only dist_calc_points and the_geom are used
            global_dist_calc_points (float): Default node interdistance.

        Returns:
            Nodes with data in the following columns:
            - id: 0-based counter generated here
            - coordinates
            - content_type: ContentType.TYPE_V2_CHANNEL
            - content_pk:  # the 0-based index into Channels (not the id!)
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

        return Nodes(
            id=np.arange(idx.size),
            coordinates=pygeos.get_coordinates(points),
            content_type=ContentType.TYPE_V2_CHANNEL,
            content_pk=idx,
        )

    def get_grid(self, nodes, channel_node_offset, connection_node_offset):
        """Compute the lines for this channel

        Args:
            self (Channels):
                connection_node_start_id and connection_node_end_id are used
            nodes (dict): additional channel nodes (see interpolate_nodes)
            channel_node_offset (int): the index of the first channel node in the
                target node array (for the lines)
            connection_node_offset (int): the index of the first connection node
                in the target node array (for the lines)

        Returns:
            Grid with data in the following columns:
            - nodes.* (see interpolate_nodes)
            - lines.id: 0-based counter generated here
            - lines.line: lines between connetion nodes and added channel
              nodes. The indices are offset using the respective parameters.
        """
        # start with the easy ones: channels that connect 2 connection nodes
        # without interpolated nodes in between
        lines_start = (
            np.vstack([self.connection_node_start_id, self.connection_node_end_id])
            + connection_node_offset
        )

        n_nodes = len(nodes)
        if n_nodes == 0:
            # if there are no interpolated nodes then we're done
            lines = Lines(id=np.arange(lines_start.shape[1]), line=lines_start.T)
            return Grid(nodes=nodes, lines=lines)

        # generate the lines that interconnect interpolated nodes
        line_ids = (
            np.vstack([np.arange(n_nodes), np.arange(1, n_nodes + 1)])
            + channel_node_offset
        )

        # connect the last line of each channel to the right
        # connection_node_end_id (instead of the next channel)
        is_channel_end = np.append(
            nodes.content_pk[1:] != nodes.content_pk[:-1], [True]
        )
        line_ids[1][is_channel_end] = (
            self.connection_node_end_id[nodes.content_pk[is_channel_end]]
            + connection_node_offset
        )

        # connect the line endings in 'lines_start' to a first interpolated
        # node (if there are interpolated nodes)
        is_channel_start = np.roll(is_channel_end, 1)
        channels_with_interp = nodes.content_pk[is_channel_start]
        lines_start[1][channels_with_interp] = (
            np.where(is_channel_start)[0] + channel_node_offset
        )


        lines = Lines(id=np.arange(lines_start.shape[1]), line=lines_start.T)
        lines += Lines(
            id=np.arange(
                lines_start.shape[1], lines_start.shape[1] + line_ids.shape[1]
            ),
            line=line_ids.T,
        )
        return Grid(nodes=nodes, lines=lines)
