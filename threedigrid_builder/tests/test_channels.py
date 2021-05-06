from numpy.testing import assert_array_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import Channels
from threedigrid_builder.grid import linear
from unittest import mock

import numpy as np
import pygeos
import pytest


@pytest.fixture
def channels():
    return Channels(
        the_geom=pygeos.linestrings([[(0, 0), (6, 0), (6, 6)]]),
        dist_calc_points=np.array([5.0]),
        id=np.array([1]),
        code=np.array(["one"]),
        connection_node_start_id=np.array([21]),
        connection_node_end_id=np.array([42]),
        calculation_type=np.array([2]),
    )


@mock.patch.object(linear.BaseLinear, "interpolate_nodes")
def test_interpolate_nodes(interpolate_nodes_m, channels):
    interpolate_nodes_m.return_value = Nodes(id=[0, 1])

    nodes = channels.interpolate_nodes(2, foo="bar")

    interpolate_nodes_m.assert_called_with(2, foo="bar")

    assert nodes is interpolate_nodes_m.return_value
    assert_array_equal(nodes.content_type, ContentType.TYPE_V2_CHANNEL)


@mock.patch.object(linear.BaseLinear, "get_lines")
def test_get_lines(get_lines_m, channels):
    get_lines_m.return_value = Lines(id=[0, 1])

    lines = channels.get_lines(2, foo="bar")

    get_lines_m.assert_called_with(2, foo="bar")

    assert lines is get_lines_m.return_value
    assert_array_equal(lines.content_type, ContentType.TYPE_V2_CHANNEL)
