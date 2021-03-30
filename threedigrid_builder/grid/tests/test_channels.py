from numpy.testing import assert_array_equal
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.grid import Channels
from threedigrid_builder.grid import ConnectionNodes

import itertools
import numpy as np
import pygeos
import pytest


@pytest.fixture
def connection_nodes():
    # Used to map connection_node_start/end_id to an index (sequence id)
    return ConnectionNodes(id=np.array([21, 25, 33, 42]))


@pytest.fixture
def one_channel():
    return Channels(
        the_geom=pygeos.linestrings([[(0, 0), (6, 0), (6, 6)]]),
        dist_calc_points=np.array([5.0]),
        id=np.array([1]),
        code=np.array(["one"]),
        connection_node_start_id=np.array([21]),
        connection_node_end_id=np.array([42]),
        calculation_type=np.array([2]),
    )


@pytest.fixture
def two_channels():
    return Channels(
        the_geom=pygeos.linestrings(
            [[(0, 0), (10, 0), (10, 10)], [(0, 0), (0, 100), (100, 100)]]
        ),
        dist_calc_points=np.array([5.0, np.nan]),
        id=np.array([1, 2]),
        code=np.array(["one", "two"]),
        connection_node_start_id=np.array([21, 25]),
        connection_node_end_id=np.array([42, 33]),
        calculation_type=np.array([2, 1]),
    )


@pytest.mark.parametrize(
    "dist,size,expected",
    [
        (3.0, 12 / 4, [(3, 0), (6, 0), (6, 3)]),  # 12 / 3  = 4 segments
        (4.0, 12 / 3, [(4, 0), (6, 2)]),  # 12 / 4  = 3 segments
        (5.0, 12 / 2, [(6, 0)]),  # 12 / 5  = 2.4 -> 2 segments
        (6.0, 12 / 2, [(6, 0)]),  # 12 / 6  = 2 segments
        (8.0, 12 / 2, [(6, 0)]),  # 12 / 8  = 1.5 -> 2
        (9.0, 12, np.empty((0, 2), dtype=float)),  # 12 / 9  = 1.33 -> 1
        (100.0, 12, np.empty((0, 2), dtype=float)),
    ],
)
def test_interpolate_nodes_one_channel(dist, size, expected, one_channel):
    one_channel.dist_calc_points[0] = dist
    nodes, segment_size = one_channel.interpolate_nodes(
        itertools.count(start=2), global_dist_calc_points=74.0
    )

    assert_array_equal(nodes.id, range(2, 2 + len(expected)))
    assert_array_equal(nodes.coordinates, expected)
    assert_array_equal(nodes.content_pk, 1)
    assert_array_equal(nodes.node_type, NodeType.NODE_1D_NO_STORAGE)
    assert_array_equal(nodes.calculation_type, 2)
    assert_array_equal(nodes.ds1d, np.arange(1, len(expected) + 1) * size)
    assert_array_equal(segment_size, size)


def test_interpolate_nodes_two_channels(two_channels):
    nodes, segment_size = two_channels.interpolate_nodes(
        itertools.count(start=2), global_dist_calc_points=50.0
    )

    expected_points = [(5, 0), (10, 0), (10, 5), (0, 50), (0, 100), (50, 100)]

    assert_array_equal(nodes.id, range(2, 8))
    assert_array_equal(nodes.coordinates, expected_points)
    assert_array_equal(nodes.content_pk, [1, 1, 1, 2, 2, 2])
    assert_array_equal(nodes.node_type, NodeType.NODE_1D_NO_STORAGE)
    assert_array_equal(nodes.calculation_type, [2, 2, 2, 1, 1, 1])
    assert_array_equal(segment_size, [5.0, 50.0])


def test_get_lines(connection_nodes, two_channels):
    nodes = Nodes(id=[10, 11, 12], content_pk=[1, 2, 2])

    lines = two_channels.get_lines(
        connection_nodes,
        nodes,
        itertools.count(start=0),
        segment_size=np.array([23, 101]),
        connection_node_offset=100,
    )

    expected_line = [(100, 10), (10, 103), (101, 11), (11, 12), (12, 102)]

    assert_array_equal(lines.id, range(5))
    assert_array_equal(lines.line, expected_line)
    assert_array_equal(lines.content_pk, [1, 1, 2, 2, 2])
    assert_array_equal(lines.content_type, ContentType.TYPE_V2_CHANNEL)
    assert_array_equal(lines.kcu, [2, 2, 1, 1, 1])
    assert_array_equal(lines.ds1d, [23, 23, 101, 101, 101])


@pytest.mark.parametrize(
    "channel_idx,expected",
    [
        ([], [(0, 3)]),
        ([1], [(0, 4), (4, 3)]),
        ([1, 1, 1], [(0, 4), (4, 5), (5, 6), (6, 3)]),
    ],
)
def test_get_lines_one_channel(channel_idx, expected, connection_nodes, one_channel):
    nodes = Nodes(id=range(4, 4 + len(channel_idx)), content_pk=channel_idx)
    lines = one_channel.get_lines(connection_nodes, nodes, itertools.count(start=0))

    assert_array_equal(lines.line, expected)


@pytest.mark.parametrize(
    "channel_idx,expected",
    [
        ([], [(0, 3), (1, 2)]),
        ([1], [(0, 4), (4, 3), (1, 2)]),
        ([1, 1, 1], [(0, 4), (4, 5), (5, 6), (6, 3), (1, 2)]),
        ([2], [(0, 3), (1, 4), (4, 2)]),
        ([2, 2, 2], [(0, 3), (1, 4), (4, 5), (5, 6), (6, 2)]),
        ([1, 2, 2], [(0, 4), (4, 3), (1, 5), (5, 6), (6, 2)]),
        ([1, 1, 2], [(0, 4), (4, 5), (5, 3), (1, 6), (6, 2)]),
    ],
)
def test_get_lines_two_channels(channel_idx, expected, connection_nodes, two_channels):
    nodes = Nodes(id=range(4, 4 + len(channel_idx)), content_pk=channel_idx)

    lines = two_channels.get_lines(connection_nodes, nodes, itertools.count(start=0))

    assert_array_equal(lines.line, expected)
