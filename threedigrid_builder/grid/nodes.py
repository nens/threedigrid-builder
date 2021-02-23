from enum import Enum, unique
from typing import Tuple
from threedi_modelchecker.threedi_model import constants
import numpy as np

from .core import array_of


@unique
class NodeType(Enum):
    NODE_2D_OPEN_WATER = 1
    NODE_2D_GROUNDWATER = 2
    NODE_1D_NO_STORAGE = 3
    NODE_1D_STORAGE = 4
    NODE_2D_BOUNDARIES = 5
    NODE_2D_GROUNDWATER_BOUNDARIES = 6
    NODE_1D_BOUNDARIES = 7


@unique
class ContentType(Enum):
    TYPE_V2_PIPE = 1
    TYPE_V2_CHANNEL = 2
    TYPE_V2_CULVERT = 3
    TYPE_V2_ORIFICE = 4
    TYPE_V2_WEIR = 5
    TYPE_V2_MANHOLE = 6
    TYPE_V2_PUMPSTATION = 7
    TYPE_V2_LEVEE = 8
    TYPE_V2_1D_LATERAL = 9
    TYPE_V2_1D_BOUNDARY_CONDITIONS = 10
    TYPE_V2_2D_BOUNDARY_CONDITIONS = 11
    TYPE_V2_CONNECTION_NODES = 12
    TYPE_V2_BREACH = 13
    TYPE_V2_CROSS_SECTION_DEF = 14
    TYPE_V2_CROSS_SECTION_LOCATION = 15
    TYPE_V2_ADDED_CALCULATION_POINT = 16
    TYPE_V2_WINDSHIELD = 17


class Node:
    id: int
    node_type: NodeType
    calculation_type: constants.CalculationType
    content_type: ContentType
    content_pk: int
    coordinates: Tuple[float, float, float, float]
    bounds: float
    bottom_level: float  # dmax
    bottom_level_groundwater: float  # dimp
    nodk: int
    nodm: int
    nodn: int
    storage_area: float


@array_of(Node)
class Nodes:
    pass
