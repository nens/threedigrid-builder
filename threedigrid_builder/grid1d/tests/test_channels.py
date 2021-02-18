import pytest
import pygeos
import numpy as np
from numpy.testing import assert_array_equal

from threedigrid_builder.grid1d import Channels


@pytest.fixture
def one_channel():
    return Channels(
        the_geom=pygeos.linestrings([[(0, 0), (6, 0), (6, 6)]]),
        dist_calc_points=np.array([5.0]),
        id=np.array([1]),
        code=np.array(["one"]),
        connection_node_start_id=np.array([21]),
        connection_node_end_id=np.array([42]),
        calculation_type=np.array([101]),
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
        calculation_type=np.array([101, 102]),
    )


@pytest.mark.parametrize(
    "dist,expected",
    [
        (3.0, pygeos.points([(3, 0), (6, 0), (6, 3)])),  # 12 / 3  = 4 segments
        (4.0, pygeos.points([(4, 0), (6, 2)])),  # 12 / 4  = 3 segments
        (5.0, pygeos.points([(6, 0)])),  # 12 / 5  = 2.4 -> 2 segments
        (6.0, pygeos.points([(6, 0)])),  # 12 / 6  = 2 segments
        (8.0, pygeos.points([(6, 0)])),  # 12 / 8  = 1.5 -> 2
        (9.0, np.empty((0,), dtype=object)),  # 12 / 9  = 1.33 -> 1
        (100.0, np.empty((0,), dtype=object)),
    ],
)
def test_interpolate_nodes(dist, expected, one_channel):
    one_channel.dist_calc_points[0] = dist
    actual = one_channel.interpolate_nodes(global_dist_calc_points=74.0)

    assert_array_equal(actual["geometry"], expected)
    assert_array_equal(actual["_channel_idx"], 0)


def test_interpolate_nodes_2(two_channels):
    actual = two_channels.interpolate_nodes(global_dist_calc_points=50.0)

    expected_points = pygeos.points(
        [(5, 0), (10, 0), (10, 5), (0, 50), (0, 100), (50, 100)]
    )
    assert_array_equal(actual["geometry"], expected_points)
    assert_array_equal(actual["_channel_idx"], [0, 0, 0, 1, 1, 1])


@pytest.mark.parametrize(
    "channel_idx,expected",
    [
        ([], [(1021, 1042)]),
        ([0], [(1021, 100), (100, 1042)]),
        ([0, 0, 0], [(1021, 100), (100, 101), (101, 102), (102, 1042)]),
    ],
)
def test_get_network(channel_idx, expected, one_channel):
    nodes = {
        "geometry": pygeos.points(np.random.random((len(channel_idx), 2))),
        "_channel_idx": np.array(channel_idx),
    }
    actual = one_channel.get_network(
        nodes, channel_node_offset=100, connection_node_offset=1000
    )
    assert_array_equal(actual, np.array(expected).T)


@pytest.mark.parametrize(
    "channel_idx,expected",
    [
        ([], [(121, 142), (125, 133)]),
        ([0], [(121, 0), (0, 142), (125, 133)]),
        ([0, 0, 0], [(121, 0), (0, 1), (1, 2), (2, 141), (125, 133)]),
        ([1], [(121, 142), (125, 0), (1, 133)]),
        ([1, 1, 1], [(121, 142), (125, 0), (0, 1), (1, 2), (2, 133)]),
        ([0, 1, 1], [(121, 0), (0, 142), (125, 1), (1, 2), (2, 133)]),
        ([0, 0, 1], [(121, 0), (0, 1), (1, 142), (125, 2), (2, 133)]),
    ],
)
def test_get_network_2(channel_idx, expected, two_channels):
    nodes = {
        "geometry": pygeos.points(np.random.random((len(channel_idx), 2))),
        "_channel_idx": np.array(channel_idx),
    }
    actual = two_channels.get_network(
        nodes, channel_node_offset=0, connection_node_offset=100
    )
    assert_array_equal(actual, np.array(expected).T)
