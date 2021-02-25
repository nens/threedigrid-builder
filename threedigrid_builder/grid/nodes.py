from threedigrid_builder.base import array_of
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import NodeType
from typing import Tuple


__all__ = ["Nodes"]


class Node:
    id: int
    node_type: NodeType
    calculation_type: CalculationType
    content_type: ContentType
    content_pk: int
    coordinates: Tuple[float, float]
    bounds: Tuple[float, float, float, float]
    dmax: float  # bottom level
    dimp: float  # bottom level groundwater
    nodk: int  # quadtree grid coordinates
    nodm: int  # quadtree grid coordinates
    nodn: int  # quadtree grid coordinates
    storage_area: float


@array_of(Node)
class Nodes:
    """Calculation node."""
