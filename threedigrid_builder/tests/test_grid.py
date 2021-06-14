from numpy.testing import assert_array_equal
from threedigrid_builder.base import GridSettings
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.base import TablesSettings
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import InitializationType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import Grid
from threedigrid_builder.grid import GridMeta
from threedigrid_builder.grid import QuadtreeStats
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


@pytest.fixture
def meta():
    return GridMeta(
        epsg_code=12432634,
        model_name="test-name",
        grid_settings=GridSettings(
            use_2d=True,
            use_1d_flow=True,
            use_2d_flow=True,
            grid_space=20.0,
            dist_calc_points=25.0,
            kmax=4,
        ),
        tables_settings=TablesSettings(
            table_step_size=0.05,
            frict_coef=0.03,
            frict_coef_type=9,
        ),
    )


@pytest.fixture
def grid2d(meta):
    quadtree_stats = QuadtreeStats(
        **{
            "lgrmin": 2,
            "kmax": 1,
            "mmax": np.array([1]),
            "nmax": np.array([1]),
            "dx": np.array([2.0]),
            "dxp": 0.5,
            "x0p": 10.0,
            "y0p": 10.0,
        }
    )
    return Grid(
        nodes=Nodes(
            id=[0, 1],
            node_type=NodeType.NODE_2D_OPEN_WATER,
            bounds=[(0, 0, 1, 1), (1, 0, 2, 1)],
        ),
        lines=Lines(id=[0]),
        meta=meta,
        quadtree_stats=quadtree_stats,
    )


@pytest.fixture
def grid1d(meta):
    return Grid(nodes=Nodes(id=[2, 3]), lines=Lines(id=[1]), meta=meta)


@pytest.mark.parametrize(
    "setting,expected_true",
    [
        ("interception_type", "has_interception"),
        ("groundwater_hydro_connectivity_type", "has_groundwater_flow"),
        ("infiltration_rate_type", "has_simple_infiltration"),
        ("groundwater_impervious_layer_level_type", "has_groundwater"),
        ("interflow_type", "has_interflow"),
    ],
)
def test_from_meta(meta, setting, expected_true):
    setattr(meta.tables_settings, setting, InitializationType.GLOBAL)
    grid = Grid.from_meta(
        epsg_code=1234,
        model_name="test",
        grid_settings=meta.grid_settings,
        tables_settings=meta.tables_settings,
    )
    for attr in (
        "has_interception",
        "has_groundwater_flow",
        "has_simple_infiltration",
        "has_groundwater",
        "has_interflow",
    ):
        assert getattr(grid.meta, attr) is (attr == expected_true)


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
    assert grid.meta == grid1d.meta
    assert grid.quadtree_stats == grid2d.quadtree_stats
    assert_array_equal(grid.nodes.id[0:2], grid2d.nodes.id)
    assert_array_equal(grid.nodes.id[2:], grid1d.nodes.id)


def test_from_channels():
    connection_nodes = mock.Mock()
    channels = mock.Mock()
    counter = mock.Mock()
    connection_node_offset = mock.Mock()
    nodes = Nodes(id=[])
    lines = Lines(id=[])

    channels.interpolate_nodes.return_value = nodes
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
        connection_node_offset=connection_node_offset,
    )


@mock.patch("threedigrid_builder.grid.cross_sections.compute_weights")
def test_set_channel_weights(compute_weights):
    # set an input grid and mock the compute_weights return value
    nodes = Nodes(
        id=[1, 2, 3],
        content_type=[ContentType.TYPE_V2_CHANNEL, -9999, ContentType.TYPE_V2_CHANNEL],
        content_pk=[1, 1, 3],
        ds1d=[2.0, 12.0, 21.0],
    )
    lines = Lines(
        id=[1, 2, 3],
        content_type=[ContentType.TYPE_V2_CHANNEL, -9999, ContentType.TYPE_V2_CHANNEL],
        content_pk=[1, 1, 3],
        ds1d=[2.0, 12.0, 21.0],
    )
    grid = Grid(nodes=nodes, lines=lines)
    compute_weights.return_value = [0, 1], [1, 2], [0.2, 0.4]  # cross1, cross2, weights
    locations = mock.Mock()
    channels = mock.Mock()

    # execute the method
    grid.set_channel_weights(locations, channels)

    # compute_weights was called correctly
    assert compute_weights.call_count == 2
    args, kwargs = compute_weights.call_args
    assert_array_equal(args[0], [1, 3])  # channel_id
    assert_array_equal(args[1], [2.0, 21.0])  # ds
    assert args[2] is locations
    assert args[3] is channels

    # node attributes cross1, cross2, cross_weight are adapted correctly
    assert_array_equal(nodes.cross1, [0, -9999, 1])
    assert_array_equal(nodes.cross2, [1, -9999, 2])
    assert_array_equal(nodes.cross_weight, [0.2, np.nan, 0.4])

    # line attributes cross1, cross2, cross_weight are adapted correctly
    assert_array_equal(lines.cross1, [0, -9999, 1])
    assert_array_equal(lines.cross2, [1, -9999, 2])
    assert_array_equal(lines.cross_weight, [0.2, np.nan, 0.4])


@mock.patch("threedigrid_builder.grid.connection_nodes.set_calculation_types")
def test_set_calculation_types(set_calculation_types, grid):
    grid.set_calculation_types()
    set_calculation_types.assert_called_with(grid.nodes, grid.lines)


@mock.patch("threedigrid_builder.grid.cross_sections.interpolate")
@mock.patch("threedigrid_builder.grid.connection_nodes.set_bottom_levels")
@mock.patch.object(Lines, "set_bottom_levels", new=mock.Mock())
@mock.patch("threedigrid_builder.grid.cross_sections.fix_dpumax")
def test_set_bottom_levels(fix_dpumax, cn_compute, cs_interpolate):
    # set an input grid and mock the compute_weights return value
    nodes = Nodes(
        id=[1, 2, 3, 4, 5],
        content_type=[
            ContentType.TYPE_V2_CHANNEL,
            -9999,
            ContentType.TYPE_V2_CHANNEL,
            ContentType.TYPE_V2_PIPE,
            ContentType.TYPE_V2_CULVERT,
        ],
        content_pk=[1, 1, 3, 2, 5],
        ds1d=[2.0, 12.0, 21.0, 15.0, 0.5],
        dmax=[np.nan, 12.0, np.nan, np.nan, np.nan],
        cross1=[5, -9999, 7, -9999, -9999],
        cross2=[6, -9999, 8, -9999, -9999],
        cross_weight=[0.2, np.nan, 0.8, np.nan, np.nan],
    )
    lines = Lines(id=[])
    cs_interpolate.return_value = [42.0, 43.0]
    grid = Grid(nodes=nodes, lines=lines)
    locations = mock.Mock()
    channels = mock.Mock()
    pipes = mock.Mock()
    pipes.compute_bottom_level.return_value = [44.0]
    weirs = mock.Mock()
    orifices = mock.Mock()
    culverts = mock.Mock()
    culverts.compute_bottom_level.return_value = [45.0]

    grid.set_bottom_levels(locations, channels, pipes, weirs, orifices, culverts)

    # cross section interpolate was called correctly
    args, _ = cs_interpolate.call_args
    assert_array_equal(args[0], [5, 7])  # cross1
    assert_array_equal(args[1], [6, 8])  # cross2
    assert_array_equal(args[2], [0.2, 0.8])  # weights
    assert args[3] is locations
    assert args[4] == "reference_level"

    # pipes.compute_bottom_level was called correctly
    args, _ = pipes.compute_bottom_level.call_args
    assert_array_equal(args[0], [2])  # pipe_id
    assert_array_equal(args[1], [15.0])  # ds

    # culverts.compute_bottom_level was called correctly
    args, _ = culverts.compute_bottom_level.call_args
    assert_array_equal(args[0], [5])  # culvert_id
    assert_array_equal(args[1], [0.5])  # ds

    # node attribute dmax is adapted correctly
    assert_array_equal(nodes.dmax, [42.0, 12.0, 43.0, 44.0, 45.0])

    # connection node set_bottom_levels was called correctly
    cn_compute.assert_called_with(
        grid.nodes, grid.lines, locations, channels, pipes, weirs, orifices, culverts
    )

    # lines set_bottom_levels was called correctly
    lines.set_bottom_levels.assert_called_with(grid.nodes, allow_nan=True)

    # fix_dpumax was called correctly
    fix_dpumax.assert_called_with(grid.lines, grid.nodes, locations)


@pytest.mark.parametrize(
    "node_coordinates,expected_lines",
    [
        ([(0.5, 0.5)], [(7, 0)]),  # first cell, center
        ([(0, 0.5)], [(7, 0)]),  # first cell, left edge
        ([(0.5, 1)], [(7, 0)]),  # first cell, top edge
        ([(0.5, 0)], [(7, 0)]),  # first cell, bottom edge
        ([(0, 1)], [(7, 0)]),  # first cell, topleft corner
        ([(0, 0)], [(7, 0)]),  # first cell, bottomleft corner
        ([(1.5, 0.5)], [(7, 1)]),  # second cell, center
        ([(2, 0.5)], [(7, 1)]),  # second cell, right edge
        ([(1.5, 1)], [(7, 1)]),  # second cell, top edge
        ([(1.5, 0)], [(7, 1)]),  # second cell, bottom edge
        ([(2, 1)], [(7, 1)]),  # second cell, topright corner
        ([(2, 0)], [(7, 1)]),  # second cell, bottomright corner
        ([(1, 1)], [(7, 0)]),  # edge between: top corner
        ([(1, 0)], [(7, 0)]),  # edge between: bottom corner
        ([(1, 0.5)], [(7, 0)]),  # edge between: middle
        ([(-1e-7, 0.5)], np.empty((0, 2), dtype=int)),  # marginally outside
        ([(2.0001, 1.5)], np.empty((0, 2), dtype=int)),  # marginally outside
        ([(1, 1.0001)], np.empty((0, 2), dtype=int)),  # marginally outside
        ([(1, -1e-7)], np.empty((0, 2), dtype=int)),  # marginally outside
        ([(0.5, 0.5), (0.5, 0.9)], [(7, 0), (8, 0)]),  # two cells, same
        ([(0.5, 0.5), (1.5, 0.5)], [(7, 0), (8, 1)]),  # two cells, different
    ],
)
def test_1d2d(node_coordinates, expected_lines, grid2d):
    grid2d.nodes += Nodes(
        id=[7, 8][: len(node_coordinates)],
        coordinates=node_coordinates,
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        calculation_type=CalculationType.CONNECTED,
    )
    grid2d.lines = Lines(id=[])

    connection_nodes = mock.Mock()
    connection_nodes.get_1d2d_properties.return_value = 0, 0
    channels = mock.Mock()
    channels.get_1d2d_properties.return_value = 0, 0
    pipes = mock.Mock()
    pipes.get_1d2d_properties.return_value = 0, 0
    locations = mock.Mock()
    culverts = mock.Mock()

    grid2d.add_1d2d(
        connection_nodes,
        channels,
        pipes,
        locations,
        culverts,
        line_id_counter=itertools.count(),
    )

    assert_array_equal(grid2d.lines.line, expected_lines)


def test_1d2d_multiple(grid2d):
    CN = ContentType.TYPE_V2_CONNECTION_NODES
    CH = ContentType.TYPE_V2_CHANNEL
    C1 = CalculationType.CONNECTED
    C2 = CalculationType.DOUBLE_CONNECTED
    grid2d.nodes += Nodes(
        id=[2, 3, 5, 7],
        coordinates=[(0.5, 0.5)] * 4,  # all the same, geo-stuff is tested elsewhere
        content_type=[CN, CN, CH, CN],
        calculation_type=[C1, C2, C2, C1],
    )
    grid2d.lines = Lines(id=[])

    connection_nodes = mock.Mock()
    connection_nodes.get_1d2d_properties.return_value = ([True, True, False], [1, 2, 3])
    channels = mock.Mock()
    channels.get_1d2d_properties.return_value = (False, [5])
    pipes = mock.Mock()
    pipes.get_1d2d_properties.return_value = 0, 0
    locations = mock.Mock()
    culverts = mock.Mock()

    grid2d.add_1d2d(
        connection_nodes,
        channels,
        pipes,
        locations,
        culverts,
        line_id_counter=itertools.count(),
    )

    args, _ = connection_nodes.get_1d2d_properties.call_args
    assert args[0] is grid2d.nodes
    assert_array_equal(args[1], [2, 3, 5])  # node_idx (offset by 2 because of 2d cells)

    args, _ = channels.get_1d2d_properties.call_args
    assert args[0] is grid2d.nodes
    assert_array_equal(args[1], [4])  # node_idx (offset by 2 because of 2d cells)
    assert args[2] is locations

    # the kcu comes from the "has_storage" from get_1d2d_properties and the calc type
    assert_array_equal(
        grid2d.lines.kcu,
        [
            LineType.LINE_1D2D_SINGLE_CONNECTED_SEWERAGE,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_SEWERAGE,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_SEWERAGE,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER,
            LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER,
        ],
    )
    # the dpumax comes from the get_1d2d_properties
    assert_array_equal(grid2d.lines.dpumax, [1.0, 2.0, 2.0, 5.0, 5.0, 3.0])
