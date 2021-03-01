from .connection_nodes import ConnectionNodes
from threedigrid_builder.base import array_of
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType, FrictionType
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


class CrossSectionLocation:
    id: int
    code: str
    the_geom: pygeos.Geometry
    definition_id: id  # refers to CrossSectionDefinition
    channel_id: id  # refers to Channel
    reference_level: float
    bank_level: float
    friction_type: FrictionType
    friction_value: float


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

        return Nodes(
            id=range(len(points)),
            coordinates=pygeos.get_coordinates(points),
            content_type=ContentType.TYPE_V2_CHANNEL,
            content_pk=self.index_to_id(idx),
        )

    def get_lines(self, connection_nodes, nodes, node_id_offsets):
        """Compute the grid (nodes + lines) for the channels.

        Fields connection_node_start_id and connection_node_end_id are used.

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
            nodes (Nodes): additional channel nodes (see interpolate_nodes)
            node_id_offsets (dict): offsets to give node indices, per type
                e.g. {ConnectionNodes: 0, Channels: 105}

        Returns:
            Lines with data in the following columns:
            - id: 0-based counter generated here
            - line: lines between connetion nodes and added channel
              nodes. The indices are offset using the respective parameters.
            - content_type: ContentType.TYPE_V2_CHANNEL
            - content_pk: the id of the Channel from which this line originates
        """
        # convert connection_node_start_id / connection_node_end_id to index
        cn_start_idx = connection_nodes.id_to_index(
            self.connection_node_start_id
        ) + node_id_offsets.get(ConnectionNodes, 0)
        cn_end_idx = connection_nodes.id_to_index(
            self.connection_node_end_id
        ) + node_id_offsets.get(ConnectionNodes, 0)

        # start with the easy ones: channels that connect 2 connection nodes
        # without interpolated nodes in between
        lines = Lines(
            id=range(len(self)),
            line=np.array([cn_start_idx, cn_end_idx]).T,
            content_pk=self.id,
            content_type=ContentType.TYPE_V2_CHANNEL,
        )

        # if there are no interpolated nodes then we're done
        if len(nodes) == 0:
            return lines

        # generate the lines that interconnect interpolated nodes
        line_ids = np.array([np.arange(len(nodes)), np.arange(1, len(nodes) + 1)])
        line_ids += node_id_offsets.get(Channels, 0)

        # connect the last line of each channel to the corresponding
        # connection_node_end_id (instead of the next channel)
        is_channel_end = np.append(
            nodes.content_pk[1:] != nodes.content_pk[:-1], [True]
        )
        channel_id = nodes.content_pk[is_channel_end]
        channel_idx = self.id_to_index(channel_id)
        line_ids[1][is_channel_end] = cn_end_idx[channel_idx]

        # connect the line endings in 'lines.line' to a first interpolated
        # node (if there are interpolated nodes)
        is_channel_start = np.roll(is_channel_end, 1)
        channel_id = nodes.content_pk[is_channel_start]
        channel_idx = self.id_to_index(channel_id)
        lines.line[channel_idx, 1] = np.where(is_channel_start)[
            0
        ] + node_id_offsets.get(Channels, 0)

        lines += Lines(
            id=range(len(lines), len(lines) + len(nodes)),
            line=line_ids.T,
            content_pk=nodes.content_pk,
            content_type=ContentType.TYPE_V2_CHANNEL,
        )
        return lines
