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
import pytest


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
        ([(0.5, 0.5), (0.5, 0.9)], [(7, 0), (8, 0)]),  # two cells, same
        ([(0.5, 0.5), (1.5, 0.5)], [(7, 0), (8, 1)]),  # two cells, different
    ],
)
def test_1d2d(node_coordinates, expected_lines, grid2d, empty_connected_points):
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
    culverts.get_1d2d_properties.return_value = 0, 0

    grid2d.add_1d2d(
        empty_connected_points,
        connection_nodes,
        channels,
        pipes,
        locations,
        culverts,
        line_id_counter=itertools.count(),
    )

    assert_array_equal(grid2d.lines.line, expected_lines)


@pytest.mark.parametrize(
    "node_coordinates", [(-1e-7, 0.5), (2.0001, 1.5), (1, 1.0001), (1, -1e-7)]
)
def test_1d2d_no_cell(node_coordinates, grid2d, empty_connected_points):
    grid2d.nodes += Nodes(
        id=[7],
        coordinates=[node_coordinates],
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        calculation_type=CalculationType.CONNECTED,
    )
    grid2d.lines = Lines(id=[])

    with pytest.raises(SchematisationError, match=".*outside of the 2D.*"):
        grid2d.add_1d2d(empty_connected_points, *((mock.Mock(),) * 6))


def test_1d2d_multiple(grid2d, empty_connected_points):
    CN = ContentType.TYPE_V2_CONNECTION_NODES
    CH = ContentType.TYPE_V2_CHANNEL
    PIPE = ContentType.TYPE_V2_PIPE
    CV = ContentType.TYPE_V2_CULVERT
    C1 = CalculationType.CONNECTED
    C2 = CalculationType.DOUBLE_CONNECTED
    grid2d.nodes += Nodes(
        id=[2, 3, 5, 7, 9, 11, 13],
        coordinates=[(0.5, 0.5)] * 7,  # all the same, geo-stuff is tested elsewhere
        content_type=[CN, CN, CH, CN, PIPE, PIPE, CV],
        calculation_type=[C1, C2, C2, C1, C1, C2, C1],
    )
    grid2d.lines = Lines(id=[])

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

    grid2d.add_1d2d(
        empty_connected_points,
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
        grid2d.lines.kcu,
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
        grid2d.lines.dpumax, [1.0, 2.0, 2.0, 5.0, 5.0, 3.0, 6.0, 7.0, 7.0, 8.0]
    )


CN = ContentType.TYPE_V2_CONNECTION_NODES
CH = ContentType.TYPE_V2_CHANNEL
PI = ContentType.TYPE_V2_PIPE
CV = ContentType.TYPE_V2_CULVERT
MH = ContentType.TYPE_V2_MANHOLE
BC = ContentType.TYPE_V2_1D_BOUNDARY_CONDITIONS


@pytest.fixture
def grid1d():
    nodes = Nodes(
        id=range(10),
        content_type=[CN, CN, CN, CN, CH, CH, CH, PI, PI, CV],
        content_pk=[1, 2, 3, 4, 1, 1, 2, 9, 9, 8],
        manhole_id=[-9999, -9999, 2, 1] + [-9999] * 6,
        boundary_id=[-9999, 1, -9999, -9999] + [-9999] * 6,
        node_type=[NodeType.NODE_1D_NO_STORAGE, NodeType.NODE_1D_BOUNDARIES]
        + [NodeType.NODE_1D_NO_STORAGE] * 8,
    )
    channels = Channels(
        id=[1, 2], connection_node_start_id=[1, 2], connection_node_end_id=[2, 1]
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