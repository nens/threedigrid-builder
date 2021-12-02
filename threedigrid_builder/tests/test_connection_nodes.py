from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.grid import Channels
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import CrossSectionLocations
from threedigrid_builder.grid import Culverts
from threedigrid_builder.grid import Orifices
from threedigrid_builder.grid import Pipes
from threedigrid_builder.grid import Weirs
from threedigrid_builder.grid.connection_nodes import set_bottom_levels
from threedigrid_builder.grid.connection_nodes import set_calculation_types

import itertools
import numpy as np
import pygeos
import pytest


@pytest.fixture
def connection_nodes():
    return ConnectionNodes(
        the_geom=pygeos.points([(0, 0), (10, 0), (10, 20), (10, 30)]),
        id=np.array([1, 3, 4, 9]),
        code=np.array(["one", "two", "", ""]),
        storage_area=np.array([15, 0, np.nan, np.nan]),
        calculation_type=np.array([1, 2, 5, 2]),
        bottom_level=np.array([2.1, np.nan, np.nan, np.nan]),
        manhole_id=[1, -9999, -9999, -9999],
        drain_level=[1.2, np.nan, np.nan, np.nan],
    )


def test_get_nodes(connection_nodes):
    counter = itertools.count(start=2)

    nodes = connection_nodes.get_nodes(counter)
    assert isinstance(nodes, Nodes)

    assert_array_equal(nodes.id, [2, 3, 4, 5])
    assert next(counter) == 6
    assert_array_equal(nodes.coordinates, [(0, 0), (10, 0), (10, 20), (10, 30)])
    assert_array_equal(nodes.content_pk, [1, 3, 4, 9])
    assert_array_equal(
        nodes.node_type,
        [
            NodeType.NODE_1D_STORAGE,
            NodeType.NODE_1D_NO_STORAGE,
            NodeType.NODE_1D_NO_STORAGE,
            NodeType.NODE_1D_NO_STORAGE,
        ],
    )
    assert_array_equal(nodes.calculation_type, connection_nodes.calculation_type)
    assert_array_equal(nodes.dmax, connection_nodes.bottom_level)
    assert_array_equal(nodes.storage_area, connection_nodes.storage_area)


@pytest.mark.parametrize(
    "kcu,expected",
    [
        ([-9999], CalculationType.ISOLATED),
        ([-9999, -9999], CalculationType.ISOLATED),
        ([LineType.LINE_1D_CONNECTED], CalculationType.CONNECTED),
        ([LineType.LINE_1D_EMBEDDED], CalculationType.EMBEDDED),
        ([LineType.LINE_1D_ISOLATED], CalculationType.ISOLATED),
        ([LineType.LINE_1D_CONNECTED, -9999], CalculationType.CONNECTED),
        ([-9999, LineType.LINE_1D_EMBEDDED], CalculationType.EMBEDDED),
        (
            [LineType.LINE_1D_ISOLATED, LineType.LINE_1D_CONNECTED],
            CalculationType.ISOLATED,
        ),
        (
            [LineType.LINE_1D_CONNECTED, LineType.LINE_1D_EMBEDDED],
            CalculationType.CONNECTED,
        ),
    ],
)
@pytest.mark.parametrize(
    "content_type",
    [
        ContentType.TYPE_V2_CHANNEL,
        ContentType.TYPE_V2_PIPE,
        ContentType.TYPE_V2_CULVERT,
    ],
)
def test_set_calculation_types_single_node(kcu, expected, content_type):
    nodes = Nodes(id=[1], content_type=ContentType.TYPE_V2_CONNECTION_NODES)
    lines = Lines(
        id=range(len(kcu)),
        content_type=content_type,
        line=[(1, 9999)] * len(kcu),
        kcu=kcu,
    )

    set_calculation_types(nodes, lines)

    assert nodes.calculation_type[0] == expected


@pytest.mark.parametrize(
    "content_type",
    [
        ContentType.TYPE_V2_CHANNEL,
        ContentType.TYPE_V2_PIPE,
        ContentType.TYPE_V2_CULVERT,
    ],
)
def test_set_calculation_types_multiple_nodes(content_type):
    nodes = Nodes(
        id=[1, 2, 3],
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        calculation_type=[-9999, -9999, CalculationType.EMBEDDED],
    )
    lines = Lines(
        id=[1, 2, 3],
        content_type=content_type,
        line=[(1, 2), (2, 9999), (9999, 1)],
        kcu=[LineType.LINE_1D_CONNECTED, -9999, LineType.LINE_1D_ISOLATED],
    )

    set_calculation_types(nodes, lines)

    assert nodes.calculation_type[0] == CalculationType.ISOLATED
    assert nodes.calculation_type[1] == CalculationType.CONNECTED
    assert nodes.calculation_type[2] == CalculationType.EMBEDDED


@pytest.mark.parametrize("_type", [Pipes, Culverts, Channels, Weirs, Orifices])
@pytest.mark.parametrize(
    "line,invert_levels,expected",
    [
        (np.empty((0, 2), dtype=int), (3.0, 4.0), np.nan),  # no line at all
        ([(2, 3)], (3.0, 4.0), np.nan),  # no line to the specific node
        ([(1, 2)], (3.0, 4.0), 3.0),  # starting point
        ([(2, 1)], (3.0, 4.0), 4.0),  # end point
        ([(1, 2), (2, 1)], (3.0, 4.0), 3.0),  # both end and start; start is lower
        ([(1, 2), (2, 1)], (4.0, 3.0), 3.0),  # both end and start; end is lower
    ],
)
def test_set_bottom_levels_single_node(line, invert_levels, expected, _type):
    nodes = Nodes(id=[1], content_type=ContentType.TYPE_V2_CONNECTION_NODES)
    lines = Lines(
        id=range(len(line)),
        content_type=_type.content_type,
        content_pk=2,
        line=line,
        invert_level_start_point=invert_levels[0],
        invert_level_end_point=invert_levels[1],
    )

    set_bottom_levels(nodes, lines)

    # assert the resulting value of dmax
    assert_almost_equal(nodes.dmax, expected)


@pytest.mark.parametrize("_type", [Pipes, Culverts, Channels, Weirs, Orifices])
def test_set_bottom_levels_multiple_nodes(_type):
    nodes = Nodes(
        id=[1, 2, 3],
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        dmax=[np.nan, np.nan, -24.0],
    )
    lines = Lines(
        id=[1, 2, 3],
        content_type=[_type.content_type, -9999, _type.content_type],
        content_pk=[2, -9999, 3],
        line=[(1, 2), (1, 2), (1, 3)],
        invert_level_start_point=[3.0, np.nan, 4.0],
        invert_level_end_point=[8.0, np.nan, 8.0],
    )

    set_bottom_levels(nodes, lines)

    # assert the resulting value of dmax
    assert_almost_equal(nodes.dmax, [3.0, 8.0, -24.0])


@pytest.mark.parametrize("structure_type", [Weirs, Orifices, Pipes, Culverts, Channels])
def test_bottom_levels_above_invert_level(structure_type, caplog):
    nodes = Nodes(
        id=[1, 2],
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        dmax=[10.0, 20.0],
        content_pk=[52, 22],
        manhole_id=[3, 6],
    )
    lines = Lines(
        id=[1],
        content_type=structure_type.content_type,
        content_pk=2,
        line=[[1, 2]],
        invert_level_start_point=3.0,
        invert_level_end_point=20.0,
    )

    # setting the bottom levels emits a warning
    set_bottom_levels(nodes, lines)
    assert caplog.messages[0].startswith("Manholes [3] have a bottom_level")

    # assert the resulting value of dmax
    assert_almost_equal(nodes.dmax, [3.0, 20.0])


def test_1d2d_properties(connection_nodes):
    """Test setup:

    - node 0: this is a manhole, drain_level is copied
    - node 1: this is skipped (not in node_idx)
    - node 2: this is not a manhole, it has 2 channels, lowest bank_level is taken
    - node 3: this is not a manhole, it has no channels: outcome is NaN
    - node 4: this is not a manhole, it 1 channel with no bank level: outcome is NaN
    """
    nodes = Nodes(
        id=[0, 1, 2, 3, 4],
        content_pk=[1, 3, 99, 4, 9],
        dmax=[1.0, 3.0, 2.0, 4.0, 2.3],
    )
    node_idx = np.array([0, 1, 3, 4])

    # channels & cs locations are so that:
    # - channel 32 (CN 1 -> 3, N 0 -> 1), start is 4.0 and end is 0.0
    # - channel 33 (CN 3 -> 9, N 1 -> 4), start&end are 1.0
    # - channel 34 (CN 3 -> 9, N 1 -> 4), 1 bank_level is nan so start&end are nan
    channels = Channels(
        id=[32, 33, 34],
        connection_node_start_id=[1, 3, 3],
        connection_node_end_id=[3, 9, 9],
        the_geom=pygeos.linestrings(
            [
                [(0, 0), (10, 0)],
                [(10, 20), (10, 30)],
                [(10, 20), (10, 30)],
            ]
        ),
    )
    locations = CrossSectionLocations(
        id=range(5),
        channel_id=[32, 32, 33, 34, 34],
        the_geom=pygeos.points([(0, 0), (10, 0), (10, 25), (10, 22), (10, 28)]),
        bank_level=[4.0, 0.0, 1.0, -10.0, np.nan],
    )

    is_closed, dpumax = connection_nodes.get_1d2d_properties(
        nodes, node_idx, channels, locations
    )

    assert_array_equal(is_closed, [True, False, False, False])
    assert_array_equal(dpumax, [1.2, 0.0, np.nan, 1.0])
