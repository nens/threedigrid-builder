from .cross_section_locations import compute_bottom_level
from threedigrid_builder.base import array_of
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType

import itertools
import logging
import numpy as np
import pygeos


logger = logging.getLogger(__name__)

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
    shape: str  # enum with classes "00", "01", "02"
    width: float
    initial_waterlevel: float
    display_name: str
    zoom_category: int


@array_of(ConnectionNode)
class ConnectionNodes:
    def get_nodes(self, node_id_counter):
        """Convert connection nodes to a nodes instance

        Args:
            node_id_counter (iterable): an iterable yielding integers

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
            - manhole_id: id of associated manhole.
            - drain_level: drain_level of associated manhole.
            - storage_area: area of connection_node.
        """
        node_type = np.where(
            self.storage_area >= 0,
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
            manhole_id=self.manhole_id,
            drain_level=self.drain_level,
            storage_area=self.storage_area,
            display_name=self.display_name,
            zoom_category=self.zoom_category,
        )
        return nodes

    def get_1d2d_properties(self, nodes, node_idx, channels, locations):
        """Compute properties (is_closed, dpumax) of 1D-2D connection node flowlines.

        Args:
            nodes (Nodes): All nodes
            node_idx (array of int): indices into nodes for which to compute properties
            channels (Channels): required for nodes without drain_level
            locations (CrossSectionLocations): to take drain_level from for nodes
                without drain_level but with channel(s)

        Returns:
            tuple of:
            - is_closed (array of bool): based on self.storage_area
            - dpumax (array of float): based on self.drain_level or locations.bank_level
        """
        # get the corresponding connection_node ids and indexes
        connection_node_id = nodes.content_pk[node_idx]
        connection_node_idx = self.id_to_index(connection_node_id)
        is_manhole = self.manhole_id[connection_node_idx] != -9999
        is_manhole_idx = connection_node_idx[is_manhole]

        # Check if manhole has drain_level below bottom_level
        has_lower_drn_lvl = (
            self.drain_level[is_manhole_idx] < self.bottom_level[is_manhole_idx]
        )
        if np.any(has_lower_drn_lvl):
            ids = self.manhole_id[is_manhole_idx[has_lower_drn_lvl]]
            logger.warning(
                f"Manholes {sorted(ids.tolist())} have a "
                f"bottom_level that is above drain_level."
            )

        # for open water nodes, compute drain level from channel bank levels
        open_water = node_idx[~is_manhole]
        dpumax = np.full(len(self), np.nan)  # easier to initialize for all conn. nodes
        for name in ("connection_node_start_id", "connection_node_end_id"):
            has_node = np.isin(getattr(channels, name), nodes.content_pk[open_water])
            cn_idx_with_channel = self.id_to_index(getattr(channels, name)[has_node])
            if name == "connection_node_start_id":
                ds = 0.0
            else:
                ds = pygeos.length(channels.the_geom[has_node])
            drain_level = compute_bottom_level(
                channels.id[has_node],
                ds,
                locations,
                channels,
                "bank_level",
            )
            _put_if_less(dpumax, cn_idx_with_channel, drain_level)
        # filter out connection nodes that were not in node_idx (non-connected ones)
        dpumax = dpumax[connection_node_idx]

        # for manholes: put in the drain level
        dpumax[is_manhole] = self.drain_level[is_manhole_idx]

        is_closed = self.storage_area[connection_node_idx] >= 0
        return is_closed, dpumax


def set_calculation_types(nodes, lines):
    """Set the calculation types for connection nodes that do not yet have one.

    The calculation_type of a connection nodes is based on
      - connection_node.manhole.calculation_type (set in ConnectionNode.get_nodes)
      - if not present, the calculation_type of adjacent channels, pipes, and culverts.
        are taken. if these are different, the precedence is:
          EMBEDDED > ISOLATED > DOUBLE_CONNECTED > CONNECTED
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
        LineType.LINE_1D_EMBEDDED,
        LineType.LINE_1D_ISOLATED,
        LineType.LINE_1D_DOUBLE_CONNECTED,
        LineType.LINE_1D_CONNECTED,
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
    # filter nan values
    is_nan = np.isnan(v)
    if is_nan.any():
        v = v[~is_nan]
        ind = ind[~is_nan]
    # same as v < a, except for how NaN is handled:
    is_less = ~(v > a[ind])
    if is_less.any():
        np.put(a, ind[is_less], v[is_less])


def set_bottom_levels(nodes, lines):
    """Set the bottom level (dmax) for connection nodes that do not yet have one.

    The bottom level (dmax) of a connection nodes is based on:
      - from the manhole.bottom_level
        - must be lower than the invert level of pipes (if not: error)
      - if not present: the lowest of all neigboring objects. we use
        lines.invert_level_start_point to get levels of channels, pipes, culverts,
        weirs, and orifices.
        - (not implemented) if no storage: invert levels should not differ more than
          cross section height! -> check this in threedi-modelchecker

    Args:
        nodes (Nodes): the nodes, those with content_type == TYPE_V2_CONNECTION_NODES
          and without a dmax will get a new dmax
        lines (Lines): the lines, including channels, pipes, weirs, and culverts
    """
    # The connection node dmax will be the lowest of these object types:
    OBJECT_TYPES = [
        ContentType.TYPE_V2_CHANNEL,
        ContentType.TYPE_V2_PIPE,
        ContentType.TYPE_V2_CULVERT,
        ContentType.TYPE_V2_WEIR,
        ContentType.TYPE_V2_ORIFICE,
    ]

    is_connection_node = nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES
    is_manhole = nodes.manhole_id != -9999
    # Copy the dmax of the manholes for later reference
    manhole_dmax = nodes.dmax[is_manhole].copy()
    # Get ids of the relevant nodes (including manholes for checking bottom levels)
    node_id = nodes.index_to_id(np.where(is_connection_node)[0])

    # line_idx: relevant lines
    line_idx = np.where(np.isin(lines.content_type, OBJECT_TYPES))[0]
    # line_idx_1: lines for which the connection node is the start
    line_idx_1 = line_idx[np.isin(lines.line[line_idx, 0], node_id)]
    if line_idx_1.size > 0:
        # get the dmax from the invert_level_start
        dmax = lines.invert_level_start_point[line_idx_1]
        # find the nodes that match to these lines and put the dmax
        _node_idx = nodes.id_to_index(lines.line[line_idx_1, 0])
        _put_if_less(nodes.dmax, _node_idx, dmax)
    # line_idx_2: lines for which the connection node is the end
    line_idx_2 = line_idx[np.isin(lines.line[line_idx, 1], node_id)]
    if line_idx_2.size > 0:
        # get the dmax from the invert_level_end
        dmax = lines.invert_level_end_point[line_idx_2]
        # find the nodes that match to these pipe lines and put the dmax
        _node_idx = nodes.id_to_index(lines.line[line_idx_2, 1])
        _put_if_less(nodes.dmax, _node_idx, dmax)

    # Check if the new node dmax is below the original manhole dmax
    has_lower_dmax = nodes.dmax[is_manhole] < manhole_dmax
    if np.any(has_lower_dmax):
        ids = nodes.manhole_id[is_manhole][has_lower_dmax]
        logger.warning(
            f"Manholes {sorted(ids.tolist())} have a "
            f"bottom_level that is above one ore more of the following connected "
            f"objects: channel reference level, pipe/culvert invert level, "
            f"weir/orifice crest level."
        )
