from enum import IntEnum
from enum import unique


@unique
class CalculationType(IntEnum):
    # channels have +100 in the calculation type, this is mapped on read
    BOUNDARY_NODE = -1
    EMBEDDED = 0
    ISOLATED = 1
    CONNECTED = 2
    BROAD_CRESTED = 3  # only orifices + weirs, corresponds
    SHORT_CRESTED = 4  # only orifices + weirs
    DOUBLE_CONNECTED = 5


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
class LineType(IntEnum):  # for kcu (calculation_type of a line)
    LINE_1D_EMBEDDED = int(CalculationType.EMBEDDED)
    LINE_1D_ISOLATED = int(CalculationType.ISOLATED)
    LINE_1D_CONNECTED = int(CalculationType.CONNECTED)
    LINE_1D_LONG_CRESTED = int(CalculationType.BROAD_CRESTED)
    LINE_1D_SHORT_CRESTED = int(CalculationType.SHORT_CRESTED)
    LINE_1D_DOUBLE_CONNECTED = int(CalculationType.DOUBLE_CONNECTED)
    LINE_1D2D_SINGLE_CONNECTED_CLOSED = 51
    LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER = 52
    LINE_1D2D_DOUBLE_CONNECTED_CLOSED = 53
    LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER = 54
    LINE_1D2D_POSSIBLE_BREACH = 55
    LINE_1D2D_ACTIVE_BREACH = 56
    LINE_1D2D_GROUNDWATER_OPEN_WATER = 57
    LINE_1D2D_GROUNDWATER_SEWER = 58  # diff to 57?
    LINE_2D_U = 98  # for internal use
    LINE_2D_V = 99  # for internal use
    LINE_2D = 100
    LINE_2D_OBSTACLE = 101  # levee
    LINE_2D_OBSTACLE_U = 102  # levee
    LINE_2D_OBSTACLE_V = 103  # levee
    LINE_2D_VERTICAL = 150
    LINE_2D_GROUNDWATER = -150
    LINE_2D_BOUNDARY_WEST = 200
    LINE_2D_BOUNDARY_EAST = 300
    LINE_2D_BOUNDARY_SOUTH = 400
    LINE_2D_BOUNDARY_NORTH = 500
    LINE_1D_BOUNDARY = 600  # for internal use


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


@unique
class FrictionType(IntEnum):
    CHEZY = 1
    MANNING = 2  # pipes have 4 here in the input, but this is mapped to 2 on read


@unique
class CrossSectionShape(IntEnum):
    CLOSED_RECTANGLE = 0  # --> only user input, convert to tabulated rectangle
    RECTANGLE = 1
    CIRCLE = 2
    EGG = 3  # --> only user input, convert to tabulated trapezium
    TABULATED_RECTANGLE = 5
    TABULATED_TRAPEZIUM = 6


@unique
class SewerageType(IntEnum):
    COMBINED = 0
    STORMWATER = 1
    WASTEWATER = 2
    TRANSPORT = 3
    OVERFLOW = 4
    SINKER = 5
    STORAGE = 6
    STORAGE_SETTLING_TANK = 7


@unique
class InitializationType(IntEnum):
    MAX = 0  # file present
    MIN = 1  # file present
    AVERAGE = 2  # file present
    NO_AGG = 3  # file present
    GLOBAL = 9  # no file present, use global value


@unique
class BoundaryType(IntEnum):
    WATERLEVEL = 1
    VELOCITY = 2
    DISCHARGE = 3
    SOMMERFELD = 5
