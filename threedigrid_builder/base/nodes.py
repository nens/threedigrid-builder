from .array import array_of
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import NodeType
from typing import Tuple

import numpy as np


__all__ = ["Nodes"]


class Node:
    id: int
    code: str
    display_name: str
    node_type: NodeType
    calculation_type: CalculationType
    content_type: ContentType
    content_pk: int
    coordinates: Tuple[float, float]
    ds1d: float  # position (arclength) along a 1D element
    bounds: Tuple[float, float, float, float]  # cell_coords in gridadmin
    pixel_coords: Tuple[int, int, int, int]  # pixel_coords in gridadmin
    dmax: float  # bottom_level or z_coordinate (?) in gridadmin
    dimp: float  # bottom level groundwater
    nodk: int  # quadtree grid coordinate z
    nodm: int  # quadtree grid coordinate x
    nodn: int  # quadtree grid coordinate y
    storage_area: float
    cross1: int  # the id of the cross section location
    cross2: int  # the id of the cross section location
    cross_weight: float


@array_of(Node)
class Nodes:
    """Calculation node."""

    def get_extent_1d(self):
        is_1d = np.isin(
            self.node_type, (NodeType.NODE_1D_NO_STORAGE, NodeType.NODE_1D_STORAGE)
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
