from numpy.testing import assert_array_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import ContentType
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


@pytest.fixture
def grid2d():
    quadtree_stats = {
        "lgrmin": 2,
        "kmax": 1,
        "mmax": np.array([1]),
        "nmax": np.array([1]),
        "dx": np.array([2.0]),
        "dxp": 0.5,
        "x0p": 10.0,
        "y0p": 10.0,
    }
    return Grid(
        nodes=Nodes(id=[0, 1]), lines=Lines(id=[0]), quadtree_stats=quadtree_stats
    )


@pytest.fixture
def grid1d():
    return Grid(nodes=Nodes(id=[2, 3]), lines=Lines(id=[1]), epsg_code=4326)


def test_from_quadtree():

    quadtree = mock.Mock()
    area_mask = mock.Mock()
    counter = mock.Mock()
    nodes = Nodes(id=[])
    lines = Lines(id=[])

    quadtree.origin = (0.0, 0.0)
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


def test_concatenate_grid(grid2d, grid1d):
    grid = grid2d + grid1d
    assert grid.epsg_code == grid1d.epsg_code
    assert grid.quadtree_stats == grid2d.quadtree_stats
    assert_array_equal(grid.nodes.id[0:2], grid2d.nodes.id)
    assert_array_equal(grid.nodes.id[2:], grid1d.nodes.id)


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


@mock.patch("threedigrid_builder.grid.cross_sections.compute_weights")
def test_set_channel_weights(compute_weights):
    # set an input grid and mock the compute_weights return value
    lines = Lines(
        id=[1, 2, 3],
        content_type=[ContentType.TYPE_V2_CHANNEL, -9999, ContentType.TYPE_V2_CHANNEL],
        content_pk=[1, 1, 3],
        ds1d=[2.0, 12.0, 21.0],
    )
    grid = Grid(nodes=Nodes(id=[]), lines=lines)
    compute_weights.return_value = [0, 1], [1, 2], [0.2, 0.4]  # cross1, cross2, weights
    locations = mock.Mock()
    channels = mock.Mock()

    # execute the method
    grid.set_channel_weights(locations, channels)

    # compute_weights was called correctly
    args, kwargs = compute_weights.call_args
    assert_array_equal(args[0], [1, 3])  # channel_id
    assert_array_equal(args[1], [2.0, 21.0])  # ds
    assert args[2] is locations
    assert args[3] is channels

    # line attributes cross1, cross2, cross_weight are adapted correctly
    assert_array_equal(lines.cross1, [0, -9999, 1])
    assert_array_equal(lines.cross2, [1, -9999, 2])
    assert_array_equal(lines.cross_weight, [0.2, np.nan, 0.4])


@mock.patch("threedigrid_builder.grid.connection_nodes.set_calculation_types")
def test_set_calculation_types(set_calculation_types, grid):
    grid.set_calculation_types()
    set_calculation_types.assert_called_with(grid.nodes, grid.lines)


@mock.patch("threedigrid_builder.grid.cross_sections.compute_bottom_level")
@mock.patch("threedigrid_builder.grid.connection_nodes.set_bottom_levels")
@mock.patch.object(Lines, "set_bottom_levels", new=mock.Mock())
@mock.patch("threedigrid_builder.grid.cross_sections.fix_dpumax")
def test_set_bottom_levels(fix_dpumax, cn_compute, cs_compute):
    # set an input grid and mock the compute_weights return value
    nodes = Nodes(
        id=[1, 2, 3],
        content_type=[ContentType.TYPE_V2_CHANNEL, -9999, ContentType.TYPE_V2_CHANNEL],
        content_pk=[1, 1, 3],
        ds1d=[2.0, 12.0, 21.0],
        dmax=[np.nan, 12.0, np.nan],
    )
    lines = Lines(id=[])
    cs_compute.return_value = [42.0, 43.0]
    grid = Grid(nodes=nodes, lines=lines)
    locations = mock.Mock()
    channels = mock.Mock()
    pipes = mock.Mock()
    weirs = mock.Mock()
    culverts = mock.Mock()

    grid.set_bottom_levels(locations, channels, pipes, weirs, culverts)

    # cross section compute_bottom_level was called correctly
    args, kwargs = cs_compute.call_args
    assert_array_equal(args[0], [1, 3])  # channel_id
    assert_array_equal(args[1], [2.0, 21.0])  # ds
    assert args[2] is locations
    assert args[3] is channels

    # node attribute dmax is adapted correctly
    assert_array_equal(nodes.dmax, [42.0, 12.0, 43.0])

    # connection node set_bottom_levels was called correctly
    cn_compute.assert_called_with(
        grid.nodes, grid.lines, locations, channels, pipes, weirs, culverts
    )

    # lines set_bottom_levels was called correctly
    lines.set_bottom_levels.assert_called_with(grid.nodes, allow_nan=False)

    # fix_dpumax was called correctly
    fix_dpumax.assert_called_with(grid.lines, grid.nodes, locations)
