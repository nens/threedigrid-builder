import itertools
import logging

import numpy as np
import shapely

from threedigrid_builder.base import Array, LineHalfs, Lines, Nodes, replace
from threedigrid_builder.constants import (
    CalculationType,
    ContentType,
    LineType,
    NodeType,
)

from .cross_section_locations import compute_bottom_level

logger = logging.getLogger(__name__)

__all__ = ["ConnectionNodes"]


class ConnectionNode:
    id: int
    the_geom: shapely.Geometry
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
    exchange_thickness: float
    hydraulic_conductivity_out: float
    hydraulic_conductivity_in: float


class ConnectionNodes(Array[ConnectionNode]):
    content_type = ContentType.TYPE_V2_CONNECTION_NODES

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

        return Nodes(
            id=itertools.islice(node_id_counter, len(self)),
            coordinates=shapely.get_coordinates(self.the_geom),
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
            manhole_indicator=self.manhole_indicator,
            has_groundwater_exchange=self.has_groundwater_exchange,
        )

    def is_closed(self, content_pk):
        """Whether object are 'closed' or 'open water' object.

        This is relevant for 1D-2D connections.
        """
        return self.storage_area[self.id_to_index(content_pk)] >= 0

    def get_1d2d_exchange_levels(self, content_pk, channels, locations, **kwargs):
        """Compute the exchange level (dpumax) for 1D-2D flowlines.

        For connection nodes 1D2D exchange levels are the drain levels. Where
        there is no drain level, the bank levels are taken from connected channels.

        Args:
            content_pk (array of int): object ids for which to compute levels
            channels (Channels): for the drain_levels
            locations (CrossSectionLocations): for the bank_levels

        Returns:
            exchange levels a.k.a. dpumax (array of float)
        """
        # get the corresponding connection_node ids and indexes
        connection_node_idx = self.id_to_index(content_pk)
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
        open_water = ~is_manhole
        dpumax = np.full(len(self), np.nan)  # easier to initialize for all conn. nodes
        for name in ("connection_node_start_id", "connection_node_end_id"):
            has_node = np.isin(getattr(channels, name), content_pk[open_water])
            cn_idx_with_channel = self.id_to_index(getattr(channels, name)[has_node])
            if name == "connection_node_start_id":
                ds = 0.0
            else:
                ds = shapely.length(channels.the_geom[has_node])
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
        return dpumax

    @property
    def has_groundwater_exchange(self):
        with np.errstate(invalid="ignore"):
            return (
                (self.storage_area > 0)
                & (self.exchange_thickness > 0)
                & np.isfinite(self.hydraulic_conductivity_out)
                & np.isfinite(self.hydraulic_conductivity_in)
            )


def set_calculation_types(nodes: Nodes, lines: Lines):
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
    PRIORITY = {
        LineType.LINE_1D_EMBEDDED: 0,
        LineType.LINE_1D_ISOLATED: 1,
        LineType.LINE_1D_DOUBLE_CONNECTED: 2,
        LineType.LINE_1D_CONNECTED: 3,
    }

    node_mask = (nodes.content_type == ContentType.TYPE_V2_CONNECTION_NODES) & (
        nodes.calculation_type == -9999
    )
    line_halfs = LineHalfs.for_connection_nodes(
        nodes,
        lines,
        line_mask=np.isin(
            lines.content_type,
            [
                ContentType.TYPE_V2_CHANNEL,
                ContentType.TYPE_V2_PIPE,
                ContentType.TYPE_V2_CULVERT,
            ],
        )
        & (lines.kcu != -9999),
        node_mask=node_mask,
    )
    line_halfs.reorder(
        np.lexsort([replace(line_halfs.kcu, PRIORITY), line_halfs.node_id])
    )

    calculation_type = line_halfs.first(line_halfs.kcu)

    # start off with ISOLATED
    nodes.calculation_type[node_mask] = CalculationType.ISOLATED
    # overwrite with the one derived from the line line_halfs
    nodes.calculation_type[
        nodes.id_to_index(line_halfs.get_reduce_id())
    ] = calculation_type


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


def set_bottom_levels(nodes: Nodes, lines: Lines):
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
    # Compute the lowest invert level for each node connection node dmax will be the lowest of these object types:
    line_halfs = LineHalfs.for_connection_nodes(
        nodes,
        lines,
        line_mask=np.isin(
            lines.content_type,
            [
                ContentType.TYPE_V2_CHANNEL,
                ContentType.TYPE_V2_PIPE,
                ContentType.TYPE_V2_CULVERT,
                ContentType.TYPE_V2_WEIR,
                ContentType.TYPE_V2_ORIFICE,
            ],
        ),
    )
    line_halfs.reorder_by("node_id")
    dmax_per_node = line_halfs.nanmin(line_halfs.invert_level)
    node_idx = nodes.id_to_index(line_halfs.get_reduce_id())

    # The new dmax is the minimum of the existing one and the one from above computation
    dmax = np.fmin(nodes.dmax[node_idx], dmax_per_node)

    # Check if the new node dmax is below the original manhole dmax
    is_manhole = nodes.manhole_id[node_idx] != -9999
    has_lower_dmax = dmax[is_manhole] < nodes.dmax[node_idx[is_manhole]]
    if np.any(has_lower_dmax):
        ids = nodes.manhole_id[node_idx[is_manhole][has_lower_dmax]]
        logger.warning(
            f"Manholes {sorted(ids.tolist())} have a "
            f"bottom_level that is above one ore more of the following connected "
            f"objects: channel reference level, pipe/culvert invert level, "
            f"weir/orifice crest level."
        )

    # Place computed dmax values
    nodes.dmax[node_idx] = dmax
