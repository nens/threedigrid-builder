from collections import OrderedDict
from typing import Tuple

import numpy as np

from threedigrid_builder.constants import (
    BoundaryType,
    CalculationType,
    ContentType,
    NodeType,
)

from .array import Array

__all__ = ["Nodes"]


PRETTY_NAMES = OrderedDict(
    [
        (ContentType.TYPE_V2_CONNECTION_NODES, "connection nodes"),
        (ContentType.TYPE_V2_CHANNEL, "channels"),
        (ContentType.TYPE_V2_PIPE, "pipes"),
        (ContentType.TYPE_V2_CULVERT, "culverts"),
    ]
)


class Node:
    id: int
    code: str
    display_name: str
    node_type: NodeType
    calculation_type: CalculationType
    content_type: ContentType
    content_pk: int
    coordinates: Tuple[float, float]
    s1d: float  # position (arclength) along a 1D element
    bounds: Tuple[float, float, float, float]  # cell_coords in gridadmin
    pixel_coords: Tuple[int, int, int, int]  # pixel_coords in gridadmin
    dmax: float  # bottom_level or z_coordinate (?) in gridadmin
    dimp: float  # bottom level groundwater
    nodk: int  # quadtree grid coordinate z
    nodm: int  # quadtree grid coordinate x
    nodn: int  # quadtree grid coordinate y
    storage_area: float
    embedded_in: int  # the id of the node in which this node is embedded
    boundary_id: int  # referring to the id of the boundary condition
    breach_ids: Tuple[int, int]  # referring to 1 or 2 ids of potential breaches
    has_dem_averaged: int  # Boolean attribute to tell if dem is averaged for node.
    boundary_type: BoundaryType
    manhole_id: int  # referring to the id of a manhole
    drain_level: float  # Drain level entered for manholes
    initial_waterlevel: float
    zoom_category: int
    # extra info, for manhole only:
    manhole_indicator: int
    surface_level: float
    shape: str  # enum with classes "00", "01", "02"
    width: float
    has_groundwater_exchange: bool  # internal


class Nodes(Array[Node]):
    """Calculation node."""

    @property
    def has_1d(self):
        return np.isin(
            self.node_type, (NodeType.NODE_1D_NO_STORAGE, NodeType.NODE_1D_STORAGE)
        ).any()

    @property
    def has_2d(self):
        return np.any(self.node_type == NodeType.NODE_2D_OPEN_WATER)

    @property
    def n_groundwater_cells(self):
        return np.count_nonzero(self.node_type == NodeType.NODE_2D_GROUNDWATER)

    def get_extent_1d(self):
        is_1d = np.isin(
            self.node_type,
            (
                NodeType.NODE_1D_NO_STORAGE,
                NodeType.NODE_1D_STORAGE,
                NodeType.NODE_1D_BOUNDARIES,
            ),
        )
        if not is_1d.any():
            return
        x, y = self.coordinates[is_1d].T
        extent = np.amin(x), np.amin(y), np.amax(x), np.amax(y)
        if any(np.isnan(val) for val in extent):
            raise ValueError("Not all 1D nodes have coordinates.")
        return extent

    def get_extent_2d(self):
        is_2d = self.node_type == NodeType.NODE_2D_OPEN_WATER
        if not is_2d.any():
            return
        x1, y1, x2, y2 = self.bounds[is_2d].T
        extent = np.amin(x1), np.amin(y1), np.amax(x2), np.amax(y2)
        if any(np.isnan(val) for val in extent):
            raise ValueError("Not all 2D nodes have coordinates.")
        return extent

    def format_message(self, where) -> str:
        """Create a message that describes the selected nodes."""
        if len(where) == 0:
            return []
        types = self.content_type[where]
        pks = self.content_pk[where]
        return ", ".join(
            [
                f"{pretty_name} {sorted(set(pks[types == content_type]))}"
                for content_type, pretty_name in PRETTY_NAMES.items()
                if np.any(types == content_type)
            ]
        )
