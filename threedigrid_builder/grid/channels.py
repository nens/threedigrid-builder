from threedigrid_builder.base import array_of
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType

import itertools
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
    def interpolate_nodes(self, node_id_counter, global_dist_calc_points):
        """Compute interpolated channel nodes

        Fields dist_calc_points and the_geom are used.

        Args:
            node_id_counter (iterable): an iterable yielding integers
            global_dist_calc_points (float): Default node interdistance.

        Returns:
            tuple of nodes (Nodes), segment_size (ndarray)
            nodes has data in the following columns:
            - id: 0-based counter generated here
            - coordinates
            - content_type: ContentType.TYPE_V2_CHANNEL
            - content_pk: the id of the Channel from which this node originates
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

        nodes = Nodes(
            id=itertools.islice(node_id_counter, len(points)),
            coordinates=pygeos.get_coordinates(points),
            content_type=ContentType.TYPE_V2_CHANNEL,
            content_pk=self.index_to_id(idx),
        )
        return nodes, segment_size

    def get_lines(
        self, connection_nodes, nodes, segment_size=None, connection_node_offset=0
    ):
        """Compute the grid (nodes + lines) for the channels.

        Fields connection_node_start_id and connection_node_end_id are used.

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
            nodes (Nodes): additional channel nodes (see interpolate_nodes)
            segment_size (ndarray of float): the segment size of each channel
              (see interpolate_nodes)
            connection_node_offset (int): offset to give connection node
              indices in the returned lines.line. Default 0.

        Returns:
            Lines with data in the following columns:
            - id: 0-based counter generated here
            - line: lines between connection nodes and added channel nodes
            - content_type: ContentType.TYPE_V2_CHANNEL
            - content_pk: the id of the Channel from which this line originates
            - ds1d: the arclength of the line (if segment_size is supplied)
            The lines are ordered by content_pk and then by position on the
            channel.
        """
        # convert connection_node_start_id / connection_node_end_id to index
        cn_start_idx = (
            connection_nodes.id_to_index(self.connection_node_start_id)
            + connection_node_offset
        )
        cn_end_idx = (
            connection_nodes.id_to_index(self.connection_node_end_id)
            + connection_node_offset
        )

        # start with the easy ones: channels that connect 2 connection nodes
        # without interpolated nodes in between
        lines = Lines(
            id=range(len(self)),
            line=np.array([cn_start_idx, cn_end_idx]).T,
            content_pk=self.id,
            content_type=ContentType.TYPE_V2_CHANNEL,
            ds1d=segment_size,
        )

        # if there are no interpolated nodes then we're done
        if len(nodes) == 0:
            return lines

        # generate the lines that interconnect interpolated nodes
        line_ids = np.array([nodes.id, np.roll(nodes.id, -1)])

        # map channel ids to channel index
        channel_idx = self.id_to_index(nodes.content_pk)

        # connect the last line of each channel to the corresponding
        # connection_node_end_id (instead of the next channel)
        is_channel_end = nodes.content_pk != np.roll(nodes.content_pk, -1)
        is_channel_end[-1] = True
        end_idx = channel_idx[is_channel_end]
        line_ids[1][is_channel_end] = cn_end_idx[end_idx]

        # connect the line endings in CN-CN lines to a first interpolated
        # node (if there are any)
        is_channel_start = np.roll(is_channel_end, 1)
        start_idx = channel_idx[is_channel_start]
        lines.line[start_idx, 1] = nodes.id[is_channel_start]

        lines += Lines(
            id=range(len(lines), len(lines) + len(nodes)),
            line=line_ids.T,
            content_pk=nodes.content_pk,
            content_type=ContentType.TYPE_V2_CHANNEL,
            ds1d=None if segment_size is None else segment_size[channel_idx],
        )

        # Reorder the lines so that they are sorted by [channel_id, position]
        id_before = lines.id
        lines.reorder_by(np.argsort(lines.content_pk))
        lines.id = id_before
        return lines

    def set_ds1d(self, segment_size, lines):
        """Set the ds1d attribute of lines from the per-channel segment size

        Args:
            segment_size (ndarray of float): the size of a segment per channel,
              as obtained by Channels().interpolate_nodes
            lines (Lines): the lines to set the ds1d on
        """
        is_channel = lines.content_type == ContentType.TYPE_V2_CHANNEL
        channel_ids = lines.content_pk[is_channel]
        channel_idx = self.id_to_index(channel_ids)
        lines.ds1d[is_channel] = segment_size[channel_idx]
