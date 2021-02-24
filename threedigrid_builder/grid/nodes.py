from enum import Enum, unique
from typing import Tuple
from threedi_modelchecker.threedi_model import constants
import numpy as np

from threedigrid_builder.base import array_of
from threedigrid_builder.constants import NodeType, ContentType, CalculationType


class Node:
    id: int
    node_type: NodeType
    calculation_type: CalculationType
    content_type: ContentType
    content_pk: int
    coordinates: Tuple[float, float]
    bounds: Tuple[float, float, float, float]
    bottom_level: float  # dmax
    bottom_level_groundwater: float  # dimp
    nodk: int
    nodm: int
    nodn: int
    storage_area: float


@array_of(Node)
class Nodes:
    pass
