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
def test_set_calculation_types_single_node(kcu, expected):
    nodes = Nodes(id=[1], content_type=ContentType.TYPE_V2_CONNECTION_NODES)
    lines = Lines(
        id=range(len(kcu)),
        content_type=ContentType.TYPE_V2_CHANNEL,
        line=[(1, 9999)] * len(kcu),
        kcu=kcu,
    )

    set_calculation_types(nodes, lines)

    assert nodes.calculation_type[0] == expected


def test_set_calculation_types_multiple_nodes():
    nodes = Nodes(
        id=[1, 2, 3],
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        calculation_type=[-9999, -9999, CalculationType.EMBEDDED],
    )
    lines = Lines(
        id=[1, 2, 3],
        content_type=ContentType.TYPE_V2_CHANNEL,
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
def test_set_bottom_levels_single_node(compute_bottom_level, line, dmax_mock, expected):
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
def test_set_bottom_levels_multiple_nodes(compute_bottom_level):
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