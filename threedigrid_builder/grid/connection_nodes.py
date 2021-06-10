from .cross_sections import compute_bottom_level
from threedigrid_builder.base import array_of
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.exceptions import SchematisationError

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
    manhole_indicator: int
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

    def get_1d2d_properties(self, nodes, node_idx):
        """Compute properties (is_sewerage, dpumax) of 1D-2D flowlines.

        Args:
            nodes (Nodes): All nodes
            node_idx (array of int): indices into nodes for which to compute properties

        Returns:
            tuple of:
            - is_sewerage (array of bool): based on self.manhole_id
            - dpumax (array of float): based on self.drain_level or nodes.dmax
        """
        # get the corresponding connection node indexes
        connection_node_id = nodes.content_pk[node_idx]
        connection_node_idx = self.id_to_index(connection_node_id)

        is_manhole = self.manhole_id[connection_node_idx] != -9999

        # assign bottom levels (for manhole: drain_level, otherwise: the node dmax)
        dpumax = np.where(
            is_manhole, self.drain_level[connection_node_idx], nodes.dmax[node_idx]
        )
        return is_manhole, dpumax


def set_calculation_types(nodes, lines):
    """Set the calculation types for connection nodes that do not yet have one.

    The calculation_type of a connection nodes is based on
      - connection_node.manhole.calculation_type (set in ConnectionNode.get_nodes)
      - if not present, the calculation_type of adjacent channels, pipes, and culverts.
        are taken. if these are different, the precedence is:
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
        | (lines.content_type == ContentType.TYPE_V2_CULVERT)
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
    """Replaces specified elements of an array with given values if  they are less."""
    # values can occur multiple times for the same index. sort descending by value so
    # that the lowest (latest) one will take precedence
    sorter = np.argsort(-v)
    v = v[sorter]
    ind = ind[sorter]
    # same as v < a, except for how NaN is handled:
    is_less = ~(v > a[ind])
    if is_less.any():
        np.put(a, ind[is_less], v[is_less])


def set_bottom_levels(
    nodes, lines, locations, channels, pipes, weirs, orifices, culverts
):
    """Set the bottom level (dmax) for connection nodes that do not yet have one.

    Lines are assumed to have the same direction as the channels and pipes.

    The bottom level (dmax) of a connection nodes is based on:
      - from the manhole.bottom_level
        - must be lower than the invert level of pipes (if not: error)
      - if not present: the lowest of all neigboring objects
        - channel: interpolate for channels (like channel node)
        - pipe and culvert: invert level start or end
           - (not implemented) if no storage: invert levels should not differ more than
             cross section height! -> check this in threedi-modelchecker
        - weir and orifice: crest level
        - (not implemented) culvert: invert level

    Args:
        nodes (Nodes): the nodes, those with content_type == TYPE_V2_CONNECTION_NODES
          and without a dmax will get a new dmax
        lines (Lines): the lines, including channels, pipes, weirs, and culverts
        locations (CrossSectionLocations): to interpolate over the channels
        channels (Channels): to interpolate over the channels
        pipes: to take invert levels
        weirs: to take crest levels
        orifices: to take crest levels
        culverts: to take invert levels
    """
    is_connection_node = nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES
    is_manhole = is_connection_node & np.isfinite(nodes.dmax)
    # Copy the dmax of the manholes for later reference
    manhole_dmax = nodes.dmax[is_manhole].copy()
    # Get ids of the relevant nodes (including manholes for checking bottom levels)
    node_id = nodes.index_to_id(np.where(is_connection_node)[0])

    # pipes & culverts
    for objs in (pipes, culverts):
        # line_idx: pipe lines
        line_idx = np.where(lines.content_type == objs.content_type)[0]
        # line_idx_1: pipe lines for which the connection node is the start
        line_idx_1 = line_idx[np.isin(lines.line[line_idx, 0], node_id)]
        if line_idx_1.size > 0:
            # get the dmax from the invert_level_start
            obj_id = lines.content_pk[line_idx_1]
            dmax = objs.invert_level_start_point[objs.id_to_index(obj_id)]
            # find the nodes that match to these pipe lines and put the dmax
            _node_idx = nodes.id_to_index(lines.line[line_idx_1, 0])
            _put_if_less(nodes.dmax, _node_idx, dmax)
        # line_idx_2: pipe lines for which the connection node is the end
        line_idx_2 = line_idx[np.isin(lines.line[line_idx, 1], node_id)]
        if line_idx_2.size > 0:
            # get the dmax from the invert_level_end
            obj_id = lines.content_pk[line_idx_2]
            dmax = objs.invert_level_end_point[objs.id_to_index(obj_id)]
            # find the nodes that match to these pipe lines and put the dmax
            _node_idx = nodes.id_to_index(lines.line[line_idx_2, 1])
            _put_if_less(nodes.dmax, _node_idx, dmax)

    # weirs & orifices
    for structures in (weirs, orifices):
        line_idx = np.where(lines.content_type == structures.content_type)[0]
        # line_idx_1: weir/orifice lines for which the connection node is the start
        line_idx_1 = line_idx[np.isin(lines.line[line_idx, 0], node_id)]
        if line_idx_1.size > 0:
            # get the dmax from the crest_level
            structure_id = lines.content_pk[line_idx_1]
            dmax = structures.crest_level[structures.id_to_index(structure_id)]
            # find the nodes that match to these weir lines and put the dmax
            _node_idx = nodes.id_to_index(lines.line[line_idx_1, 0])
            _put_if_less(nodes.dmax, _node_idx, dmax)
        # line_idx_2: weir/orifice lines for which the connection node is the end
        line_idx_2 = line_idx[np.isin(lines.line[line_idx, 1], node_id)]
        if line_idx_2.size > 0:
            # get the dmax from the crest_level
            structure_id = lines.content_pk[line_idx_1]
            dmax = structures.crest_level[structures.id_to_index(structure_id)]
            # find the nodes that match to these weir lines and put the dmax
            _node_idx = nodes.id_to_index(lines.line[line_idx_2, 1])
            _put_if_less(nodes.dmax, _node_idx, dmax)

    # Check if the new node dmax is below the original manhole dmax
    has_lower_dmax = nodes.dmax[is_manhole] < manhole_dmax
    if np.any(has_lower_dmax):
        ids = nodes.content_pk[is_manhole][has_lower_dmax]
        raise SchematisationError(
            f"Connection nodes {sorted(ids.tolist())} have a manhole with a "
            f"bottom_level that is above a connected pipe invert level, weir crest "
            f"level, or orifice crest level."
        )

    # Get the indices and ids of the relevant nodes (now excluding manholes)
    node_id = nodes.index_to_id(np.where(is_connection_node & ~is_manhole)[0])

    # line_idx: channel lines
    line_idx = np.where(lines.content_type == ContentType.TYPE_V2_CHANNEL)[0]
    # line_idx_1: channel lines for which the connection node is the start
    line_idx_1 = line_idx[np.isin(lines.line[line_idx, 0], node_id)]
    if line_idx_1.size > 0:
        # get dmax by interpolating along the line's channel at ds=0 (start)
        channel_id = lines.content_pk[line_idx_1]
        channel_ds = np.zeros(len(channel_id))
        dmax = compute_bottom_level(channel_id, channel_ds, locations, channels)
        # find the nodes that match to these channel lines and put the dmax
        _node_idx = nodes.id_to_index(lines.line[line_idx_1, 0])
        _put_if_less(nodes.dmax, _node_idx, dmax)
    # line_idx_2: channel lines for which the connection node is the end
    line_idx_2 = line_idx[np.isin(lines.line[line_idx, 1], node_id)]
    if line_idx_2.size > 0:
        # get dmax by interpolating along the line's channel at ds=length(channel)
        channel_id = lines.content_pk[line_idx_2]
        channel_ds = pygeos.length(channels.the_geom[channels.id_to_index(channel_id)])
        dmax = compute_bottom_level(channel_id, channel_ds, locations, channels)
        # find the nodes that match to these channel lines and put the dmax
        _node_idx = nodes.id_to_index(lines.line[line_idx_2, 1])
        _put_if_less(nodes.dmax, _node_idx, dmax)
