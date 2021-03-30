from .cross_sections import compute_dmax
from threedigrid_builder.base import array_of
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import ManholeIndicator
from threedigrid_builder.constants import NodeType

import itertools
import numpy as np
import pygeos


__all__ = ["ConnectionNodes"]


class ConnectionNode:
    id: int
    the_geom: pygeos.Geometry
    code: str
    storage_area: float
    manhole_id: int
    calculation_type: CalculationType
    manhole_indicator: ManholeIndicator
    bottom_level: float
    drain_level: float
    surface_level: float
    manhole_shape: str  # enum with classes "00", "01", "02"
    manhole_width: float


@array_of(ConnectionNode)
class ConnectionNodes:
    def get_nodes(self, node_id_counter):
        """Convert connection nodes to a nodes instance

        Args:
            node_id_counter (iterable): an iterable yielding integers
            global_dist_calc_points (float): Default node interdistance.

        Returns:
            tuple of nodes (Nodes), segment_size (ndarray)
            nodes has data in the following columns:
            - id: 0-based counter generated here
            - coordinates
            - content_type: TYPE_V2_CHANNEL
            - content_pk: the id of the Channel from which this node originates
            - node_type: NODE_1D_STORAGE or NODE_1D_NO_STORAGE depending on storage_area
            - calculation_type: from calculation_type (which comes from manhole)
            - dmax: from bottom_level (which comes from manhole)
        """
        node_type = np.where(
            self.storage_area > 0,
            int(NodeType.NODE_1D_STORAGE),
            int(NodeType.NODE_1D_NO_STORAGE),
        )
        nodes = Nodes(
            id=itertools.islice(node_id_counter, len(self)),
            coordinates=pygeos.get_coordinates(self.the_geom),
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            content_pk=self.id,
            node_type=node_type,
            calculation_type=self.calculation_type,
            dmax=self.bottom_level,
        )
        return nodes


def set_calculation_types(nodes, lines):
    """Set the calculation types for connection nodes that do not yet have one.

    The calculation_type of a connection nodes is based on
      - connection_node.manhole.calculation_type (set in ConnectionNode.get_nodes)
      - if not present, the calculation_type of adjacent channels or pipes are taken
        if these are different, the precedence is:
          ISOLATED > DOUBLE_CONNECTED > CONNECTED > EMBEDDED
      - if not present, then the calculation type becomes ISOLATED

    Args:
        nodes (Nodes): the nodes, those with content_type == TYPE_V2_CONNECTION_NODES
          and without a calculation_type will get a new calculation_type
        lines (Lines): the lines, including channels and pipes
    """
    # Get the indices of the relevant nodes and lines
    node_idx = np.where(
        (nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES)
        & (nodes.calculation_type == -9999)
    )[0]
    line_idx = np.where(
        (lines.content_type == ContentType.TYPE_V2_CHANNEL)
        | (lines.content_type == ContentType.TYPE_V2_PIPE)
    )[0]

    for _type in (
        LineType.LINE_1D_ISOLATED,
        LineType.LINE_1D_DOUBLE_CONNECTED,
        LineType.LINE_1D_CONNECTED,
        LineType.LINE_1D_EMBEDDED,
    ):
        # Do the nodes have a line with this _type?
        has_this_type = np.isin(
            nodes.index_to_id(node_idx),
            lines.line[line_idx[lines.kcu[line_idx] == _type]],
        )
        # set the type (note: LineType and CalculationType have equal enum values)
        nodes.calculation_type[node_idx[has_this_type]] = _type
        # these nodes are 'done', skip them in the next loop
        node_idx = node_idx[~has_this_type]

    # Remaining nodes get ISOLATED
    nodes.calculation_type[node_idx] = CalculationType.ISOLATED


def _put_if_less(a, ind, v):
    """Replaces specified elements of an array with given values if they are less."""
    is_less = v < a[ind]
    if is_less.any():
        np.put(a, ind[is_less], v[is_less])


def set_bottom_level(nodes, lines, locations, channels, pipes, weirs, culverts):
    """Set the bottom level (dmax) for connection nodes that do not yet have one.

    Lines are assumed to have the same direction as the channels and pipes.

    The bottom level (dmax) of a connection nodes is based on:
      - from the manhole.bottom_level
      - if not present: the lowest of all neigboring objects
        - channel: interpolate for channels (like channel node)
        - pipe: invert level (if no storage & invert levels differ by more
          than cross section height: error)
        - weir: crest level
        - culvert: invert level

    Args:
        nodes (Nodes): the nodes, those with content_type == TYPE_V2_CONNECTION_NODES
          and without a dmax will get a new dmax
        lines (Lines): the lines, including channels, pipes, weirs, and culverts
        locations (CrossSectionLocations): to interpolate over the channels
        channels (Channels): to interpolate over the channels
        pipes: to take invert levels, not yet implemented
        weirs: to take crest levels, not yet implemented
        culverts: to take invert levels, not yet implemented
    """
    # Get the indices and ids of the relevant nodes
    node_idx = np.where(
        (nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES)
        & (np.isnan(nodes.dmax))
    )[0]
    node_id = nodes.index_to_id(node_idx)

    # Connection node that are channel start / endpoints:
    line_mask = lines.content_type == ContentType.TYPE_V2_CHANNEL
    ch_start = node_idx[np.isin(node_id, lines.line[line_mask, 0])]
    ch_end = node_idx[np.isin(node_id, lines.line[line_mask, 1])]
    if ch_start.size > 0:
        channel_ids = nodes.content_pk[ch_start]
        channel_ds = 0
        dmax = compute_dmax(channel_ids, channel_ds, locations, channels)
        _put_if_less(nodes.dmax, ch_start, dmax)
    if ch_end.size > 0:
        channel_ids = nodes.content_pk[ch_end]
        channel_ds = pygeos.length(channels.the_geom[channels.id_to_index(channel_ids)])
        dmax = compute_dmax(channel_ids, channel_ds, locations, channels)
        _put_if_less(nodes.dmax, ch_end, dmax)

    # pipes
    line_mask = lines.content_type == ContentType.TYPE_V2_PIPE
    pipe_start = node_idx[np.isin(node_id, lines.line[line_mask, 0])]
    pipe_end = node_idx[np.isin(node_id, lines.line[line_mask, 1])]
    if pipe_start.size > 0:
        raise NotImplementedError()
    if pipe_end.size > 0:
        raise NotImplementedError()

    # weirs
    line_mask = lines.content_type == ContentType.TYPE_V2_WEIR
    has_weir = node_idx[np.isin(node_id, lines.line[line_mask, :])]
    if has_weir.size > 0:
        raise NotImplementedError()

    # culverts
    line_mask = lines.content_type == ContentType.TYPE_V2_CULVERT
    has_culvert = node_idx[np.isin(node_id, lines.line[line_mask, :])]
    if has_culvert.size > 0:
        raise NotImplementedError()
