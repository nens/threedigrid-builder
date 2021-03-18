from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import Grid
from unittest import mock

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


def test_from_quadtree():

    quadtree = mock.Mock()
    area_mask = mock.Mock()
    counter = mock.Mock()
    nodes = Nodes(id=[])
    lines = Lines(id=[])

    quadtree.get_nodes_lines.return_value = nodes, lines

    grid = Grid.from_quadtree(
        quadtree=quadtree,
        area_mask=area_mask,
        node_id_counter=counter,
        line_id_counter=counter,
    )

    assert isinstance(grid, Grid)
    assert grid.nodes is nodes
    assert grid.lines is lines


def test_from_connection_nodes():
    connection_nodes = mock.Mock()
    counter = mock.Mock()
    nodes = Nodes(id=[])

    connection_nodes.get_nodes.return_value = nodes

    grid = Grid.from_connection_nodes(connection_nodes, counter)

    connection_nodes.get_nodes.assert_called_with(counter)

    assert grid.nodes is nodes
    assert len(grid.lines) == 0


def test_from_channels():
    connection_nodes = mock.Mock()
    channels = mock.Mock()
    counter = mock.Mock()
    segment_size = mock.Mock()
    connection_node_offset = mock.Mock()
    nodes = Nodes(id=[])
    lines = Lines(id=[])

    channels.interpolate_nodes.return_value = nodes, segment_size
    channels.get_lines.return_value = lines
    grid = Grid.from_channels(
        connection_nodes,
        channels,
        global_dist_calc_points=100.0,
        node_id_counter=counter,
        line_id_counter=counter,
        connection_node_offset=connection_node_offset,
    )

    assert isinstance(grid, Grid)
    assert grid.nodes is nodes
    assert grid.lines is lines

    channels.interpolate_nodes.assert_called_with(counter, 100.0)
    channels.get_lines.assert_called_with(
        connection_nodes,
        nodes,
        counter,
        segment_size=segment_size,
        connection_node_offset=connection_node_offset,
    )


@mock.patch("threedigrid_builder.grid.grid.compute_weights")
def test_set_channel_weights(compute_weights, grid):
    locations = mock.Mock()
    channels = mock.Mock()

    grid.set_channel_weights(locations, channels)

    compute_weights.assert_called_with(grid.lines, locations, channels)


@mock.patch("threedigrid_builder.grid.grid.set_calculation_types")
def test_set_calculation_types(set_calculation_types, grid):
    grid.set_calculation_types()
    set_calculation_types.assert_called_with(grid.nodes, grid.lines)
