from . import geo_utils
from threedigrid_builder.base import array_of
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import NodeType

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
            - node_type: NodeType.NODE_1D_NO_STORAGE
        """
        # insert default dist_calc_points where necessary
        dists = self.dist_calc_points.copy()  # copy because of inplace edits
        dists[~np.isfinite(dists)] = global_dist_calc_points
        dists[dists <= 0] = global_dist_calc_points

        # interpolate the node geometries
        points, index = geo_utils.segmentize(self.the_geom, dists)

        # construct the nodes with available attributes
        nodes = Nodes(
            id=itertools.islice(node_id_counter, len(points)),
            coordinates=pygeos.get_coordinates(points),
            content_type=ContentType.TYPE_V2_CHANNEL,
            content_pk=self.index_to_id(index),
            node_type=NodeType.NODE_1D_NO_STORAGE,
            calculation_type=self.calculation_type[index],
        )
        return nodes

    def get_lines(
        self,
        connection_nodes,
        nodes,
        line_id_counter,
        connection_node_offset=0,
    ):
        """Compute the grid lines for the channels.

        Fields connection_node_start_id and connection_node_end_id are used.

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
            nodes (Nodes): additional channel nodes (see interpolate_nodes)
            line_id_counter (iterable): an iterable yielding integers
            connection_node_offset (int): offset to give connection node
              indices in the returned lines.line. Default 0.

        Returns:
            Lines with data in the following columns:
            - id: 0-based counter generated here
            - line: lines between connection nodes and added channel nodes
            - content_type: ContentType.TYPE_V2_CHANNEL
            - content_pk: the id of the Channel from which this line originates
            - ds1d: the arclength of the line (if segment_size is supplied)
            - kcu: the calculation_type of the Channel
            The lines are ordered by content_pk and then by position on the
            channel.
        """
        # count the number of segments per channel
        (n_lines,) = self.the_geom.shape
        node_line_idx = self.id_to_index(nodes.content_pk)
        segment_counts = np.bincount(node_line_idx, minlength=n_lines) + 1

        # cut the channel geometries into segment geometries
        start_s, end_s, segment_idx = geo_utils.segment_start_end(
            self.the_geom, segment_counts
        )
        segments = geo_utils.line_substring(self.the_geom, start_s, end_s, segment_idx)

        # set the right node indices for each segment
        first_idx, last_idx = geo_utils.counts_to_ranges(segment_counts)
        last_idx -= 1  # convert slice end into last index
        line = np.full((len(segments), 2), -9999, dtype=np.int32)

        # convert connection_node_start_id to index and put it at first segments' start
        line[first_idx, 0] = (
            connection_nodes.id_to_index(self.connection_node_start_id)
            + connection_node_offset
        )
        # convert connection_node_end_id to index and put it at last segments' end
        line[last_idx, 1] = (
            connection_nodes.id_to_index(self.connection_node_end_id)
            + connection_node_offset
        )
        # set node indices to line start where segment start is not a conn. node
        mask = np.ones(len(segments), dtype=bool)
        mask[first_idx] = False
        line[mask, 0] = nodes.id
        # set node indices to line end where segment end is not a conn. node
        mask = np.ones(len(segments), dtype=bool)
        mask[last_idx] = False
        line[mask, 1] = nodes.id

        # construct the result
        return Lines(
            id=itertools.islice(line_id_counter, len(segments)),
            line_geometries=segments,
            line=line,
            content_pk=self.id[segment_idx],
            content_type=ContentType.TYPE_V2_CHANNEL,
            ds1d=end_s - start_s,
            kcu=self.calculation_type[segment_idx],
        )
