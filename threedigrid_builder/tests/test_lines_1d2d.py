from numpy.testing import assert_array_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import Channels
from threedigrid_builder.grid import ConnectedPoints
from threedigrid_builder.grid import Culverts
from threedigrid_builder.grid import Grid
from threedigrid_builder.grid import Pipes
from unittest import mock

import itertools
import numpy as np
import pygeos
import pytest


CN = ContentType.TYPE_V2_CONNECTION_NODES
CH = ContentType.TYPE_V2_CHANNEL
PI = ContentType.TYPE_V2_PIPE
CV = ContentType.TYPE_V2_CULVERT
MH = ContentType.TYPE_V2_MANHOLE
BC = ContentType.TYPE_V2_1D_BOUNDARY_CONDITIONS
NODE_1D = NodeType.NODE_1D_NO_STORAGE
NODE_1DBC = NodeType.NODE_1D_BOUNDARIES
ISO = CalculationType.ISOLATED
C1 = CalculationType.CONNECTED
C2 = CalculationType.DOUBLE_CONNECTED


@pytest.fixture
def grid2d():
    return Grid(
        nodes=Nodes(
            id=[0, 1],
            node_type=NodeType.NODE_2D_OPEN_WATER,
            bounds=[(0, 0, 1, 1), (1, 0, 2, 1)],
        ),
        lines=Lines(id=[0]),
    )


@pytest.fixture
def empty_connected_points():
    return ConnectedPoints(id=[])


@pytest.mark.parametrize(
    "node_coordinates,expected_lines",
    [
        ([(0.5, 0.5)], [(0, 7)]),  # first cell, center
        ([(0, 0.5)], [(0, 7)]),  # first cell, left edge
        ([(0.5, 1)], [(0, 7)]),  # first cell, top edge
        ([(0.5, 0)], [(0, 7)]),  # first cell, bottom edge
        ([(0, 1)], [(0, 7)]),  # first cell, topleft corner
        ([(0, 0)], [(0, 7)]),  # first cell, bottomleft corner
        ([(1.5, 0.5)], [(1, 7)]),  # second cell, center
        ([(2, 0.5)], [(1, 7)]),  # second cell, right edge
        ([(1.5, 1)], [(1, 7)]),  # second cell, top edge
        ([(1.5, 0)], [(1, 7)]),  # second cell, bottom edge
        ([(2, 1)], [(1, 7)]),  # second cell, topright corner
        ([(2, 0)], [(1, 7)]),  # second cell, bottomright corner
        ([(1, 1)], [(0, 7)]),  # edge between: top corner
        ([(1, 0)], [(0, 7)]),  # edge between: bottom corner
        ([(1, 0.5)], [(0, 7)]),  # edge between: middle
        ([(0.5, 0.5), (0.5, 0.9)], [(0, 7), (0, 8)]),  # two cells, same
        ([(0.5, 0.5), (1.5, 0.5)], [(0, 7), (1, 8)]),  # two cells, different
    ],
)
def test_get_lines(node_coordinates, expected_lines, grid2d, empty_connected_points):
    grid2d.nodes += Nodes(
        id=[7, 8][: len(node_coordinates)],
        coordinates=node_coordinates,
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        calculation_type=CalculationType.CONNECTED,
    )

    connection_nodes = mock.Mock()
    connection_nodes.get_1d2d_properties.return_value = 0, 0
    channels = mock.Mock()
    channels.get_1d2d_properties.return_value = 0, 0
    pipes = mock.Mock()
    pipes.get_1d2d_properties.return_value = 0, 0
    locations = mock.Mock()
    culverts = mock.Mock()
    culverts.get_1d2d_properties.return_value = 0, 0

    actual_lines = empty_connected_points.get_lines(
        grid2d.cell_tree,
        grid2d.nodes,
        connection_nodes,
        channels,
        pipes,
        locations,
        culverts,
        line_id_counter=itertools.count(),
    )

    assert_array_equal(actual_lines.line, expected_lines)


@pytest.mark.parametrize(
    "node_coordinates", [(-1e-7, 0.5), (2.0001, 1.5), (1, 1.0001), (1, -1e-7)]
)
def test_get_lines_no_cell(node_coordinates, grid2d, empty_connected_points):
    grid2d.nodes += Nodes(
        id=[7],
        coordinates=[node_coordinates],
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        calculation_type=CalculationType.CONNECTED,
    )

    with pytest.raises(SchematisationError, match=".*outside of the 2D.*"):
        empty_connected_points.get_lines(
            grid2d.cell_tree, grid2d.nodes, *((mock.Mock(),) * 6)
        )


def test_get_lines_multiple(grid2d, empty_connected_points):
    grid2d.nodes += Nodes(
        id=[2, 3, 5, 7, 9, 11, 13],
        coordinates=[(0.5, 0.5)] * 7,  # all the same, geo-stuff is tested elsewhere
        content_type=[CN, CN, CH, CN, PI, PI, CV],
        calculation_type=[C1, C2, C2, C1, C1, C2, C1],
    )

    connection_nodes = mock.Mock()
    connection_nodes.get_1d2d_properties.return_value = (
        [True, True, True, False],
        [1, 2, 2, 3],
    )
    channels = mock.Mock()
    channels.get_1d2d_properties.return_value = ([False, False], [5, 5])
    pipes = mock.Mock()
    pipes.get_1d2d_properties.return_value = ([True, True, True], [6, 7, 7])
    locations = mock.Mock()
    culverts = mock.Mock()
    culverts.get_1d2d_properties.return_value = ([True], [8])

    actual_lines = empty_connected_points.get_lines(
        grid2d.cell_tree,
        grid2d.nodes,
        connection_nodes,
        channels,
        pipes,
        locations,
        culverts,
        line_id_counter=itertools.count(),
    )

    args, _ = connection_nodes.get_1d2d_properties.call_args
    assert args[0] is grid2d.nodes
    assert_array_equal(
        args[1], [2, 3, 3, 5]
    )  # node_idx (offset by 2 because of 2d cells)

    args, _ = channels.get_1d2d_properties.call_args
    assert args[0] is grid2d.nodes
    assert_array_equal(args[1], [4, 4])  # node_idx (offset by 2 because of 2d cells)
    assert args[2] is locations

    args, _ = pipes.get_1d2d_properties.call_args
    assert args[0] is grid2d.nodes
    assert_array_equal(args[1], [6, 7, 7])  # node_idx (offset by 2 because of 2d cells)
    assert args[2] is connection_nodes

    args, _ = culverts.get_1d2d_properties.call_args
    assert args[0] is grid2d.nodes
    assert_array_equal(args[1], [8])  # node_idx (offset by 2 because of 2d cells)
    assert args[2] is connection_nodes

    # the kcu comes from the "has_storage" from get_1d2d_properties and the calc type
    assert_array_equal(
        actual_lines.kcu,
        [
            LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_CLOSED,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_CLOSED,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER,
            LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER,
            LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_CLOSED,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_CLOSED,
            LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED,
        ],
    )
    # the dpumax comes from the get_1d2d_properties
    assert_array_equal(
        actual_lines.dpumax, [1.0, 2.0, 2.0, 5.0, 5.0, 3.0, 6.0, 7.0, 7.0, 8.0]
    )


def test_get_lines_multiple_with_conn_points(grid2d):
    grid2d.nodes += Nodes(
        id=[2, 3, 5, 7, 9, 11, 13],
        coordinates=[(0.5, 0.5)] * 7,  # all the same, geo-stuff is tested elsewhere
        content_type=[CN, CN, CH, CN, PI, PI, CV],
        calculation_type=[C1, C2, C2, C1, C1, C2, C1],
    )
    connected_points = ConnectedPoints(
        id=range(3),
        the_geom=pygeos.points([(1.5, 0.5), (1.5, 0.5), (1.5, 0.5)]),
        exchange_level=[np.nan, 0.1, 0.2],
        # don't set the content_pk etc., instead, we patch .get_node_index()
    )

    connection_nodes = mock.Mock()
    connection_nodes.get_1d2d_properties.return_value = True, [1, 2, 2, 3]
    channels = mock.Mock()
    channels.get_1d2d_properties.return_value = (False, [5, 5])
    pipes = mock.Mock()
    pipes.get_1d2d_properties.return_value = (True, [6, 7, 7])
    locations = mock.Mock()
    culverts = mock.Mock()
    culverts.get_1d2d_properties.return_value = (True, [8])

    with mock.patch.object(connected_points, "get_node_index") as get_node_index:
        get_node_index.return_value = np.array([2, 2, 6]) + 2
        actual_lines = connected_points.get_lines(
            grid2d.cell_tree,
            grid2d.nodes,
            connection_nodes,
            channels,
            pipes,
            locations,
            culverts,
            line_id_counter=itertools.count(),
        )

    # the dpumax comes from the get_1d2d_properties, but is overridded by exchange_level
    assert_array_equal(
        actual_lines.dpumax, [1.0, 2.0, 2.0, 5.0, 0.1, 3.0, 6.0, 7.0, 7.0, 0.2]
    )
    # the cell id is different for the 3 user-supplied connected points
    assert_array_equal(actual_lines.line[:, 0], [0, 0, 0, 1, 1, 0, 0, 0, 0, 1])


@pytest.fixture
def grid1d():
    nodes = Nodes(
        id=range(10),
        content_type=[CN, CN, CN, CN, CH, CH, CH, PI, PI, CV],
        content_pk=[1, 2, 3, 4, 1, 1, 2, 9, 9, 8],
        manhole_id=[-9999, -9999, 2, 1] + [-9999] * 6,
        boundary_id=[-9999, 1, -9999, -9999] + [-9999] * 6,
        node_type=[NODE_1D, NODE_1DBC] + [NODE_1D] * 8,
        calculation_type=[C1, C2, ISO, ISO, C1, C1, ISO, C2, C2, ISO],
    )
    channels = Channels(
        id=[1, 2, 3],
        connection_node_start_id=[1, 2, 1],
        connection_node_end_id=[2, 1, 4],
    )
    pipes = Pipes(id=[9], connection_node_start_id=[3], connection_node_end_id=[4])
    culverts = Culverts(
        id=[8], connection_node_start_id=[4], connection_node_end_id=[3]
    )
    return nodes, channels, pipes, culverts


@pytest.mark.parametrize(
    "content_type,content_pk,node_number,expected",
    [
        ([BC], [1], [-9999], [1]),
        ([MH], [2], [-9999], [2]),
        ([MH], [1], [-9999], [3]),
        ([CH], [1], [1], [0]),
        ([CH], [1], [2], [4]),
        ([CH], [1], [3], [5]),
        ([CH], [1], [4], [1]),
        ([CH], [2], [1], [1]),
        ([CH], [2], [2], [6]),
        ([CH], [2], [3], [0]),
        ([CH], [3], [1], [0]),
        ([CH], [3], [2], [3]),
        ([PI], [9], [1], [2]),
        ([PI], [9], [2], [7]),
        ([PI], [9], [3], [8]),
        ([PI], [9], [4], [3]),
        ([CV], [8], [1], [3]),
        ([CV], [8], [2], [9]),
        ([CV], [8], [3], [2]),
        ([CH] * 6, [1] * 6, [1, 2, 2, 3, 3, 4], [0, 4, 4, 5, 5, 1]),
        ([CH, PI], [1, 9], [2, 2], [4, 7]),
        ([MH, MH], [2, 1], [-9999, -9999], [2, 3]),
    ],
)
def test_get_node_index(content_type, content_pk, node_number, expected, grid1d):
    connected_points = ConnectedPoints(
        id=range(len(expected)),
        content_type=content_type,
        content_pk=content_pk,
        node_number=node_number,
    )

    actual = connected_points.get_node_index(*grid1d)

    assert_array_equal(actual, expected)


@pytest.mark.parametrize(
    "content_type,content_pk,node_number,expected",
    [
        ([BC], [2], [-9999], r".*\[10\] refer to non-existing 1D boundary conditions."),
        ([MH], [3], [-9999], r".*\[10\] refer to non-existing manholes."),
        ([CH], [4], [1], r".*\[10\] refer to non-existing channels."),
        ([PI], [3], [1], r".*\[10\] refer to non-existing pipes."),
        ([CV], [3], [1], r".*\[10\] refer to non-existing culverts."),
        ([CH], [4], [2], r".*\[10\] refer to non-existing channels."),
        ([PI], [3], [2], r".*\[10\] refer to non-existing pipes."),
        ([CV], [3], [2], r".*\[10\] refer to non-existing culverts."),
        ([CH], [1], [0], r".*\[10\] have node numbers below 1."),
        ([CH], [1], [5], r".*\[10\] have too large node numbers."),
        ([CH], [3], [3], r".*\[10\] have too large node numbers."),
    ],
)
def test_get_node_index_err(content_type, content_pk, node_number, expected, grid1d):
    connected_points = ConnectedPoints(
        id=range(len(content_type)),
        content_type=content_type,
        content_pk=content_pk,
        node_number=node_number,
        calculation_point_id=range(10, 10 + len(content_type)),
    )

    with pytest.raises(SchematisationError, match=expected):
        connected_points.get_node_index(*grid1d)


@pytest.mark.parametrize(
    "cp_node_idx,expected_line_cp_idx",
    [
        ([0, 4, 5], [0, -9999, -9999, 1, 2, -9999, -9999, -9999, -9999]),
        ([1, 5, 7], [-9999, 0, -9999, -9999, 1, 2, -9999, -9999, -9999]),
        ([1, 1, 0], [2, 0, 1, -9999, -9999, -9999, -9999, -9999, -9999]),
        ([7, 8, 7], [-9999, -9999, -9999, -9999, -9999, 0, 2, 1, -9999]),
    ],
)
def test_get_line_mappings(cp_node_idx, expected_line_cp_idx, grid1d):
    connected_points = ConnectedPoints(
        id=range(len(cp_node_idx)),
    )

    line_node_idx, line_cp_idx, line_is_double = connected_points.get_line_mappings(
        grid1d[0], cp_node_idx
    )

    assert_array_equal(line_cp_idx, expected_line_cp_idx)
    assert_array_equal(line_node_idx, [0, 1, 1, 4, 5, 7, 7, 8, 8])
    assert_array_equal(line_is_double, [False, True, True, False, False] + [True] * 4)


@pytest.mark.parametrize(
    "cp_node_idx,expected_msg",
    [
        ([2], ".*refer to objects that do not have a 'connected'.*"),
        ([0, 0], ".*are a second reference to objects that do not have a 'double.*"),
        ([1, 1, 1], ".*have too many connected points."),
    ],
)
def test_get_line_mappings_err(cp_node_idx, expected_msg, grid1d):
    connected_points = ConnectedPoints(
        id=range(len(cp_node_idx)),
        calculation_point_id=range(10, 10 + len(cp_node_idx)),
    )

    with pytest.raises(SchematisationError, match=expected_msg):
        connected_points.get_line_mappings(grid1d[0], cp_node_idx)
