from .array import array_of
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import NodeType
from typing import Tuple


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


@array_of(Node)
class Nodes:
    """Calculation node."""
