from enum import IntEnum
from enum import unique


@unique
class CalculationType(IntEnum):
    # -1 tot 99 ::= Node | 'verbinding-zonder-subnodes(="vertices")'
    BOUNDARY_NODE = -1
    EMBEDDED_NODE = 0
    STANDALONE_NODE = 1
    CONNECTED_NODE = 2
    BROAD_CRESTED = 3  # only orifices + weirs, corresponds
    SHORT_CRESTED = 4  # only orifices + weirs
    DOUBLE_CONNECTED = 5
    # 100 tot oneindig ::= 'verbinding-met-subnodes(="vertices")'
    # a.k.a. pinpoint
    EMBEDDED_VERTEX = 100
    STANDALONE_VERTEX = 101
    CONNECTED_VERTEX = 102
    DOUBLE_CONNECTED_VERTEX = 105


@unique
class NodeType(IntEnum):
    NODE_2D_OPEN_WATER = 1
    NODE_2D_GROUNDWATER = 2
    NODE_1D_NO_STORAGE = 3
    NODE_1D_STORAGE = 4
    NODE_2D_BOUNDARIES = 5
    NODE_2D_GROUNDWATER_BOUNDARIES = 6
    NODE_1D_BOUNDARIES = 7


@unique
class ContentType(IntEnum):
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
