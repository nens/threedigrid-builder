from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import Channels
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import Pipes
from threedigrid_builder.grid.connection_nodes import get_1d2d_lines
from threedigrid_builder.grid.connection_nodes import set_bottom_levels
from threedigrid_builder.grid.connection_nodes import set_calculation_types
from unittest import mock

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
        storage_area=np.array([-2, 0, 15, np.nan]),
        calculation_type=np.array([1, 2, 5, 2]),
        bottom_level=np.array([2.1, 5.2, np.nan, np.nan]),
        manhole_id=[1, 2, -9999, -9999],
        drain_level=[1.2, 3.4, np.nan, np.nan],
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
            NodeType.NODE_1D_NO_STORAGE,
            NodeType.NODE_1D_NO_STORAGE,
            NodeType.NODE_1D_STORAGE,
            NodeType.NODE_1D_NO_STORAGE,
        ],
    )
    assert_array_equal(nodes.calculation_type, connection_nodes.calculation_type)
    assert_array_equal(nodes.dmax, connection_nodes.bottom_level)


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
    "content_type", [ContentType.TYPE_V2_CHANNEL, ContentType.TYPE_V2_PIPE]
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
    "content_type", [ContentType.TYPE_V2_CHANNEL, ContentType.TYPE_V2_PIPE]
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


@pytest.mark.parametrize(
    "line,dmax_mock,expected",
    [
        (np.empty((0, 2), dtype=int), (), np.nan),  # no line at all
        ([(2, 3)], (), np.nan),  # no line to the specific node
        ([(1, 2)], (3.0,), 3.0),  # starting point
        ([(2, 1)], (3.0,), 3.0),  # end point
        ([(1, 2), (2, 1)], (3.0, 4.0), 3.0),  # both end and start; start is lower
        ([(1, 2), (2, 1)], (4.0, 3.0), 3.0),  # both end and start; end is lower
        ([(1, 2), (1, 3)], ([3.0, 4.0],), 3.0),  # two lines; first is lower
        ([(1, 2), (1, 3)], ([4.0, 3.0],), 3.0),  # two lines; last is lower
    ],
)
@mock.patch("threedigrid_builder.grid.connection_nodes.compute_bottom_level")
def test_set_bottom_levels_single_node_channels(
    compute_bottom_level, line, dmax_mock, expected
):
    nodes = Nodes(id=[1], content_type=ContentType.TYPE_V2_CONNECTION_NODES)
    lines = Lines(
        id=range(len(line)),
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=2,
        line=line,
    )
    channels = Channels(
        id=[2],
        the_geom=[pygeos.linestrings([[0, 0], [0, 10]])],
    )
    locations = mock.Mock()
    pipes = mock.Mock()
    culverts = mock.Mock()
    weirs = mock.Mock()

    compute_bottom_level.side_effect = [np.atleast_1d(x) for x in dmax_mock]
    set_bottom_levels(nodes, lines, locations, channels, pipes, weirs, culverts)

    # assert the correct call to compute_bottom_level
    assert compute_bottom_level.call_count == len(dmax_mock)

    # assert the resulting value of dmax
    assert_almost_equal(nodes.dmax, expected)


@mock.patch("threedigrid_builder.grid.connection_nodes.compute_bottom_level")
def test_set_bottom_levels_multiple_nodes_channels(compute_bottom_level):
    nodes = Nodes(
        id=[1, 2, 3],
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        dmax=[np.nan, np.nan, 24.0],
    )
    lines = Lines(
        id=[1, 2, 3],
        content_type=[ContentType.TYPE_V2_CHANNEL, -9999, ContentType.TYPE_V2_CHANNEL],
        content_pk=[2, -9999, 3],
        line=[(1, 2), (1, 2), (1, 3)],
    )
    channels = Channels(
        id=[2, 3],
        the_geom=[
            pygeos.linestrings([[0, 0], [0, 10]]),
            pygeos.linestrings([[0, 0], [0, 2]]),
        ],
    )
    locations = mock.Mock()
    pipes = mock.Mock()
    culverts = mock.Mock()
    weirs = mock.Mock()

    compute_bottom_level.side_effect = (np.array([3.0, 4.0]), np.array([8.0]))
    set_bottom_levels(nodes, lines, locations, channels, pipes, weirs, culverts)

    # assert the correct call to compute_dmax
    assert compute_bottom_level.call_count == 2
    (first_call, _), (second_call, _) = compute_bottom_level.call_args_list
    assert_array_equal(first_call[0], [2, 3])  # channel ids for channel starts
    assert_array_equal(first_call[1], [0.0, 0.0])  # ds for channel starts
    assert first_call[2] is locations
    assert first_call[3] is channels
    assert_array_equal(second_call[0], [2])  # channel ids for channel endings
    assert_array_equal(second_call[1], [10.0])  # ds for channel endings (= length)
    assert second_call[2] is locations
    assert second_call[3] is channels

    # assert the resulting value of dmax
    assert_almost_equal(nodes.dmax, [3.0, 8.0, 24.0])


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
def test_set_bottom_levels_single_node_pipes(line, invert_levels, expected):
    nodes = Nodes(id=[1], content_type=ContentType.TYPE_V2_CONNECTION_NODES)
    lines = Lines(
        id=range(len(line)),
        content_type=ContentType.TYPE_V2_PIPE,
        content_pk=2,
        line=line,
    )
    pipes = Pipes(
        id=[2],
        invert_level_start_point=[invert_levels[0]],
        invert_level_end_point=[invert_levels[1]],
    )
    locations = mock.Mock()
    channels = mock.Mock()
    culverts = mock.Mock()
    weirs = mock.Mock()

    set_bottom_levels(nodes, lines, locations, channels, pipes, weirs, culverts)

    # assert the resulting value of dmax
    assert_almost_equal(nodes.dmax, expected)


def test_set_bottom_levels_multiple_nodes_pipes():
    nodes = Nodes(
        id=[1, 2, 3],
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        dmax=[np.nan, np.nan, -24.0],
    )
    lines = Lines(
        id=[1, 2, 3],
        content_type=[ContentType.TYPE_V2_PIPE, -9999, ContentType.TYPE_V2_PIPE],
        content_pk=[2, -9999, 3],
        line=[(1, 2), (1, 2), (1, 3)],
    )
    pipes = Pipes(
        id=[2, 3],
        invert_level_start_point=[3.0, 4.0],
        invert_level_end_point=[8.0, 8.0],
    )
    locations = mock.Mock()
    channels = mock.Mock()
    culverts = mock.Mock()
    weirs = mock.Mock()

    set_bottom_levels(nodes, lines, locations, channels, pipes, weirs, culverts)

    # assert the resulting value of dmax
    assert_almost_equal(nodes.dmax, [3.0, 8.0, -24.0])


def test_bottom_levels_above_invert_level():
    nodes = Nodes(
        id=[1, 2],
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        dmax=[10.0, 20.0],
        content_pk=[52, 22],
    )
    lines = Lines(
        id=[1],
        content_type=ContentType.TYPE_V2_PIPE,
        content_pk=2,
        line=[[1, 2]],
    )
    pipes = Pipes(
        id=[2],
        invert_level_start_point=[3.0],
        invert_level_end_point=[4.0],
    )
    locations = mock.Mock()
    channels = mock.Mock()
    culverts = mock.Mock()
    weirs = mock.Mock()

    # assert the resulting value of dmax
    with pytest.raises(SchematisationError) as e:
        set_bottom_levels(nodes, lines, locations, channels, pipes, weirs, culverts)
        assert str(e).startswith("Connection nodes [22, 52] have a manhole")


@pytest.fixture
def cell_args():
    cells = pygeos.box(0, 0, 1, 1), pygeos.box(1, 0, 2, 1)
    cell_ids = np.array([5, 6])
    cell_tree = pygeos.STRtree(cells)
    return cell_tree, cell_ids


@pytest.mark.parametrize(
    "node_coordinates,expected_lines",
    [
        ([(0.5, 0.5)], [(0, 5)]),  # first cell, center
        ([(0, 0.5)], [(0, 5)]),  # first cell, left edge
        ([(0.5, 1)], [(0, 5)]),  # first cell, top edge
        ([(0.5, 0)], [(0, 5)]),  # first cell, bottom edge
        ([(0, 1)], [(0, 5)]),  # first cell, topleft corner
        ([(0, 0)], [(0, 5)]),  # first cell, bottomleft corner
        ([(1.5, 0.5)], [(0, 6)]),  # second cell, center
        ([(2, 0.5)], [(0, 6)]),  # second cell, right edge
        ([(1.5, 1)], [(0, 6)]),  # second cell, top edge
        ([(1.5, 0)], [(0, 6)]),  # second cell, bottom edge
        ([(2, 1)], [(0, 6)]),  # second cell, topright corner
        ([(2, 0)], [(0, 6)]),  # second cell, bottomright corner
        ([(1, 1)], [(0, 5)]),  # edge between: top corner
        ([(1, 0)], [(0, 5)]),  # edge between: bottom corner
        ([(1, 0.5)], [(0, 5)]),  # edge between: middle
        ([(-1e-7, 0.5)], np.empty((0, 2), dtype=int)),  # marginally outside
        ([(2.0001, 1.5)], np.empty((0, 2), dtype=int)),  # marginally outside
        ([(1, 1.0001)], np.empty((0, 2), dtype=int)),  # marginally outside
        ([(1, -1e-7)], np.empty((0, 2), dtype=int)),  # marginally outside
        ([(0.5, 0.5), (0.5, 0.9)], [(0, 5), (1, 5)]),  # two cells, same
        ([(0.5, 0.5), (1.5, 0.5)], [(0, 5), (1, 6)]),  # two cells, different
    ],
)
def test_1d2d(node_coordinates, expected_lines, connection_nodes, cell_args):
    nodes = Nodes(
        id=range(len(node_coordinates)),
        coordinates=node_coordinates,
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        calculation_type=CalculationType.CONNECTED,
    )
    actual = get_1d2d_lines(nodes, *cell_args, connection_nodes, itertools.count())
    assert_array_equal(actual.line, expected_lines)


def test_1d2d_multiple(connection_nodes, cell_args):
    CN = ContentType.TYPE_V2_CONNECTION_NODES
    CH = ContentType.TYPE_V2_CHANNEL
    C1 = CalculationType.CONNECTED
    C2 = CalculationType.DOUBLE_CONNECTED
    nodes = Nodes(
        id=[0, 2, 5, 7],
        coordinates=[(0.5, 0.5)] * 4,  # all the same, geo-stuff is tested elsewhere
        content_type=[CN, CN, CH, CN],
        content_pk=[1, 3, 1, 4],
        calculation_type=[C1, C2, C2, C1],
        dmax=[1.0, 3.0, 2.0, 4.0],
    )

    actual = get_1d2d_lines(nodes, *cell_args, connection_nodes, itertools.count())
    # all connect to the same cell, but the double-connected exists twice
    assert_array_equal(actual.line, [(0, 5), (2, 5), (2, 5), (7, 5)])
    # only the manholes (conn. nodes 1 and 3) have storage
    assert_array_equal(
        actual.kcu,
        [
            LineType.LINE_1D2D_SINGLE_CONNECTED_WITH_STORAGE,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_WITH_STORAGE,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_WITH_STORAGE,
            LineType.LINE_1D2D_SINGLE_CONNECTED_WITHOUT_STORAGE,
        ],
    )
    # for the manholes, conn_node.drain_level is copied, otherwise, dmax is taken
    assert_array_equal(actual.dpumax, [1.2, 3.4, 3.4, 4.0])
