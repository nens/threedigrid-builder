from numpy.testing import assert_array_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import Grid
from unittest import mock

import itertools
import numpy as np
import pygeos
import pytest


@pytest.fixture
def connection_nodes():
    return ConnectionNodes(
        the_geom=pygeos.points([(0, 0), (10, 0)]),
        id=np.array([1, 3]),
        code=np.array(["one", "two"]),
    )


@pytest.fixture
def grid():
    return Grid(nodes=Nodes(id=[]), lines=Lines(id=[]))


def test_from_connection_nodes(connection_nodes):
    counter = itertools.count(start=2)

    grid = Grid.from_connection_nodes(connection_nodes, counter)

    assert_array_equal(grid.nodes.id, [2, 3])
    assert next(counter) == 4
    assert_array_equal(grid.nodes.coordinates, [(0, 0), (10, 0)])
    assert_array_equal(grid.nodes.content_pk, [1, 3])


def test_from_channels():
    connection_nodes = mock.Mock()
    channels = mock.Mock()
    counter = mock.Mock()
    segment_size = mock.Mock()
    nodes = Nodes(id=[])
    lines = Lines(id=[])

    channels.interpolate_nodes.return_value = nodes, segment_size
    channels.get_lines.return_value = lines
    grid = Grid.from_channels(
        connection_nodes,
        channels,
        global_dist_calc_points=100.0,
        node_id_counter=counter,
    )

    assert isinstance(grid, Grid)
    assert grid.nodes is nodes
    assert grid.lines is lines

    channels.interpolate_nodes.assert_called_with(counter, 100.0)
    channels.get_lines.assert_called_with(
        connection_nodes, nodes, segment_size=segment_size
    )


@mock.patch("threedigrid_builder.grid.grid.cross_sections.compute_weights")
def test_set_channel_weights(compute_weights, grid):
    locations = mock.Mock()
    channels = mock.Mock()

    grid.set_channel_weights(locations, channels)

    compute_weights.assert_called_with(grid.lines, locations, channels)
