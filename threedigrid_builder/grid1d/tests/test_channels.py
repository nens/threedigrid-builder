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


@pytest.fixture
def two_channels():
    return Channels(
        the_geom=pygeos.linestrings(
            [[(0, 0), (10, 0), (10, 10)], [(0, 0), (0, 100), (100, 100)]]
        ),
        dist_calc_points=np.array([5.0, np.nan]),
        id=np.array([1, 2]),
        code=np.array(["one", "two"]),
        connection_node_start_id=np.array([21, 33]),
        connection_node_end_id=np.array([21, 42]),
        calculation_type=np.array([101, 102]),
    )


def test_interpolate_nodes_2(two_channels):
    actual = two_channels.interpolate_nodes(global_dist_calc_points=50.0)

    expected_points = pygeos.points(
        [(5, 0), (10, 0), (10, 5), (0, 50), (0, 100), (50, 100)]
    )
    assert_array_equal(actual["geometry"], expected_points)
    assert_array_equal(actual["calculation_type"][:3], 101)
    assert_array_equal(actual["calculation_type"][3:], 102)
