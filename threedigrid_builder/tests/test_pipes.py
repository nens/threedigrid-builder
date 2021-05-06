from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import Channels
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import linear
from threedigrid_builder.grid import Pipes
from threedigrid_builder.grid.pipes import compute_bottom_level
from unittest import mock

import numpy as np
import pygeos
import pytest


@pytest.fixture
def connection_nodes():
    # Used to map connection_node_start/end_id to an index (sequence id)
    return ConnectionNodes(
        id=[21, 25, 33, 42],
        the_geom=pygeos.points([(0, 21), (1, 25), (2, 33), (3, 42)]),
    )


@pytest.fixture
def pipes():
    return Pipes(
        id=[1, 2],
        dist_calc_points=[5.0, np.nan],
        code=["one", "two"],
        connection_node_start_id=[21, 21],
        connection_node_end_id=[25, 42],
        calculation_type=[2, 1],
        invert_level_start_point=[3.0, 5.0],
        invert_level_end_point=[4.0, 6.0],
    )


@pytest.fixture
def pipes_with_geom(pipes, connection_nodes):
    pipes.set_geometries(connection_nodes)
    return pipes


def test_set_geometries(pipes, connection_nodes):
    pipes.set_geometries(connection_nodes)

    expected_geometries = pygeos.linestrings(
        [
            [(0, 21), (1, 25)],
            [(0, 21), (3, 42)],
        ]
    )
    assert pygeos.equals(pipes.the_geom, expected_geometries).all()


def test_interpolate_nodes_no_geometries(pipes):
    with pytest.raises(ValueError, match=".*Call set_geometries first.*"):
        pipes.interpolate_nodes(2, foo="bar")


@mock.patch.object(linear.BaseLinear, "interpolate_nodes")
def test_interpolate_nodes(interpolate_nodes_m, pipes_with_geom):
    interpolate_nodes_m.return_value = Nodes(id=[0, 1])

    nodes = pipes_with_geom.interpolate_nodes(2, foo="bar")

    interpolate_nodes_m.assert_called_with(2, foo="bar")

    assert nodes is interpolate_nodes_m.return_value
    assert_array_equal(nodes.content_type, ContentType.TYPE_V2_PIPE)


@mock.patch.object(linear.BaseLinear, "get_lines")
def test_get_lines(get_lines_m, pipes):
    get_lines_m.return_value = Lines(id=[0, 1])

    lines = pipes.get_lines(2, foo="bar")

    get_lines_m.assert_called_with(2, foo="bar")

    assert lines is get_lines_m.return_value
    assert_array_equal(lines.content_type, ContentType.TYPE_V2_PIPE)


@pytest.mark.parametrize(
    "pipe_ids,ds,expected",
    [
        ([1], [5.0], [3.5]),
        ([2], [0.5], [5.5]),
        ([1, 2], [5.0, 0.5], [3.5, 5.5]),
        ([2, 1], [0.5, 5.0], [5.5, 3.5]),
        ([1, 1, 2, 2], [5.0, 7.5, 0.5, 0.25], [3.5, 3.75, 5.5, 5.25]),
    ],
)
def test_compute_bottom_level(pipe_ids, ds, pipes, expected):
    # set geometries with lengths 10 and 1 (resp. id 1 and 2)
    # invert levels are [3, 4] for id=1 and [5, 6] for id=2
    pipes.the_geom = pygeos.linestrings([[(0, 0), (0, 10)], [(2, 2), (3, 2)]])

    actual = compute_bottom_level(pipe_ids, ds, pipes)

    assert_almost_equal(actual, expected)


def test_compute_bottom_level_raises_no_geom(pipes):
    with pytest.raises(ValueError, match=".*Call set_geometries first.*"):
        compute_bottom_level([1], [5.0], pipes)


@pytest.mark.parametrize(
    "pipe_ids,ds",
    [
        ([1], [10.1]),
        ([2], [-1e-7]),
    ],
)
def test_compute_bottom_level_raises_out_of_bounds(pipe_ids, ds, pipes):
    pipes.the_geom = pygeos.linestrings([[(0, 0), (0, 10)], [(2, 2), (3, 2)]])
    with pytest.raises(ValueError, match=".*outside of the pipe bounds.*"):
        compute_bottom_level(pipe_ids, ds, pipes)
