from .grid import Grid  # TODO this import should not be here
from threedigrid_builder.base import array_of
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType

import numpy as np
import pygeos


__all__ = ["Channels"]


class Channel:
    id: int
    code: str
    the_geom: pygeos.Geometry
    dist_calc_points: float
    connection_node_start_id: int
    connection_node_end_id: int
    calculation_type: CalculationType


@array_of(Channel)
class Channels:
    def interpolate_nodes(self, global_dist_calc_points):
        """Compute interpolated channel nodes

        Fields dist_calc_points and the_geom are used.

        Args:
            global_dist_calc_points (float): Default node interdistance.

        Returns:
            Nodes with data in the following columns:
            - id: 0-based counter generated here
            - coordinates
            - content_type: ContentType.TYPE_V2_CHANNEL
            - content_pk: the 0-based index into Channels (not the channel id)
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
        """Compute the grid (nodes + lines) for the channels.

        Fields connection_node_start_id and connection_node_end_id are used.

        Args:
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
            - lines.content_type: ContentType.TYPE_V2_CHANNEL
            - content_pk: the 0-based index into Channels (not the channel id)
        """
        # start with the easy ones: channels that connect 2 connection nodes
        # without interpolated nodes in between
        lines_start = np.vstack(
            [
                self.connection_node_start_id,
                self.connection_node_end_id,
            ]
        )
        lines_start += connection_node_offset
        lines = Lines(
            id=np.arange(len(self)),
            line=lines_start.T,
            content_pk=np.arange(len(self)),  # indices into self (channels)
            content_type=ContentType.TYPE_V2_CHANNEL,
        )

        n_nodes = len(nodes)
        if n_nodes == 0:
            # if there are no interpolated nodes then we're done
            return Grid(nodes=nodes, lines=lines)

        # generate the lines that interconnect interpolated nodes
        line_ids = np.vstack([np.arange(n_nodes), np.arange(1, n_nodes + 1)])
        line_ids += channel_node_offset

        # connect the last line of each channel to the corresponding
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
        lines.line[channels_with_interp, 1] = (
            np.where(is_channel_start)[0] + channel_node_offset
        )

        lines += Lines(
            id=np.arange(len(lines), len(lines) + len(nodes)),
            line=line_ids.T,
            content_pk=nodes.content_pk,
            content_type=ContentType.TYPE_V2_CHANNEL,
        )
        return Grid(nodes=nodes, lines=lines)
