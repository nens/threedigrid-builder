from itertools import count
from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import Channels
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import embed_channels
from threedigrid_builder.grid import embed_nodes
from threedigrid_builder.grid import Grid
from unittest import mock

import numpy as np
import pygeos
import pytest


NODE_2D = NodeType.NODE_2D_OPEN_WATER
NODE_1D = NodeType.NODE_1D_STORAGE
EMBEDDED = CalculationType.EMBEDDED
ISOLATED = CalculationType.ISOLATED


@pytest.fixture
def grid2d():
    return Grid(
        nodes=Nodes(
            id=[0, 1, 2, 3],
            node_type=NODE_2D,
            bounds=[(0, 0, 1, 1), (1, 0, 2, 1), (0, 1, 1, 2), (1, 1, 2, 2)],
        ),
        lines=Lines(id=[0], line=[[0, 1]]),
    )


def test_embed_nodes(grid2d):
    grid = grid2d + Grid(
        nodes=Nodes(
            id=[4, 5, 6, 7],
            node_type=NODE_1D,
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            calculation_type=[EMBEDDED, ISOLATED, EMBEDDED, EMBEDDED],
            coordinates=[(1.0, 0.2), (0.9, 0.9), (1.3, 0.3), (1.6, 0.7)],
        ),
        lines=Lines(id=[1, 2, 3], line=[(4, 5), (4, 6), (5, 6)]),
    )
    embedded = embed_nodes(grid, count(10))

    assert_array_equal(grid.nodes.id, [0, 1, 2, 3, 4])
    assert_array_equal(grid.nodes.node_type, [NODE_2D] * 4 + [NODE_1D])
    assert_array_equal(grid.lines.line, [(0, 1), (0, 4), (0, 1), (4, 1)])

    assert_array_equal(embedded.id, [10, 11, 12])  # new ids
    assert_array_equal(embedded.calculation_type, EMBEDDED)
    assert_array_equal(embedded.embedded_in, [0, 1, 1])


def test_embed_node_outside_2D(grid2d):
    grid = grid2d + Grid(
        nodes=Nodes(
            id=[4, 5],
            node_type=NODE_1D,
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            content_pk=[15, 16],
            calculation_type=EMBEDDED,
            coordinates=[(0.5, 2.2), (0.5, 0.5)],
        ),
        lines=Lines(id=[]),
    )

    with pytest.raises(SchematisationError, match=r".*\[15\] are outside the 2D cells"):
        embed_nodes(grid, count(10))


def test_embed_node_two_interconnected(grid2d):
    grid = grid2d + Grid(
        nodes=Nodes(
            id=[4, 5],
            node_type=NODE_1D,
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            content_pk=[4, 5],
            calculation_type=EMBEDDED,
            coordinates=[(0.8, 0.8), (0.5, 0.5)],
        ),
        lines=Lines(id=[2], line=[(4, 5)]),
    )

    with pytest.raises(SchematisationError, match=r".*\[4, 5\] connect to.*"):
        embed_nodes(grid, count(10))


@pytest.fixture
def connection_nodes():
    # Used to map connection_node_start/end_id to an index (sequence id)
    return ConnectionNodes(
        id=[21, 25, 33, 42],
        the_geom=pygeos.points([(0, 21), (1, 25), (2, 33), (3, 42)]),
    )


def test_embed_channels_multiple(grid2d, connection_nodes):
    channels = Channels(
        id=[0, 1],
        calculation_type=EMBEDDED,
        the_geom=[
            pygeos.linestrings([(0.1, 0.1), (1.5, 0.1)]),
            pygeos.linestrings([(0.7, 0.2), (1.8, 0.2), (1.8, 0.9), (1.8, 1.8)]),
        ],
        connection_node_start_id=[21, 25],
        connection_node_end_id=[42, 33],
    )

    nodes, lines = embed_channels(
        channels,
        connection_nodes,
        grid2d.cell_tree,
        count(2),
        count(3),
        connection_node_offset=10,
    )

    assert_array_equal(nodes.id, [2])
    assert_array_equal(nodes.content_type, ContentType.TYPE_V2_CHANNEL)
    assert_array_equal(nodes.content_pk, [1])
    assert_array_equal(nodes.coordinates, [(1.8, 0.2)])
    assert_array_equal(nodes.calculation_type, EMBEDDED)
    assert_array_equal(nodes.s1d, [1.1])  # halfway the 2 velocity points (see below)
    assert_array_equal(nodes.embedded_in, [1])

    assert_array_equal(lines.id, [3, 4, 5])
    assert_array_equal(lines.line, [(10, 13), (11, 1), (1, 12)])  # connect to cell 1
    assert_array_equal(lines.content_type, ContentType.TYPE_V2_CHANNEL)
    assert_array_equal(lines.content_pk, [0, 1, 1])
    assert_array_equal(lines.kcu, EMBEDDED)
    assert_array_equal(lines.cross1, -9999)
    assert_array_equal(lines.cross2, -9999)
    assert_array_equal(lines.cross_weight, np.nan)
    assert_almost_equal(lines.s1d, [0.9, 0.3, 1.9])
    assert_almost_equal(lines.ds1d, [1.4, 1.1, 1.6])
    assert_almost_equal(
        pygeos.get_coordinates(lines.line_geometries[0]), [(0.1, 0.1), (1.5, 0.1)]
    )
    assert_almost_equal(
        pygeos.get_coordinates(lines.line_geometries[1]), [(0.7, 0.2), (1.8, 0.2)]
    )
    assert_almost_equal(
        pygeos.get_coordinates(lines.line_geometries[2]),
        [(1.8, 0.2), (1.8, 0.9), (1.8, 1.8)],
    )


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize(
    "channel,lines_s1d,embedded_in",
    [
        # ([(0.5, 0.5), (0.9, 0.5)], [], []),  # no cell crossing, errors
        ([(0.5, 0.5), (1.8, 0.5)], [0.5], []),  # horizontal, 1 crossing
        ([(0.5, 0.5), (0.5, 1.8)], [0.5], []),  # vertical, 1 crossing
        ([(0.5, 0.5), (0.5, 0.9), (1.5, 0.9)], [0.9], []),  # 1 crossing, more coords
        ([(0.5, 0.5), (1.8, 0.5), (1.8, 0.9)], [0.5], []),  # 1 crossing, more coords
        ([(0.5, 1.5), (0.5, 0.1), (1.8, 0.1)], [0.5, 1.9], [0]),  # corner, bottomleft
        ([(0.5, 0.5), (1.8, 0.5), (1.8, 1.3)], [0.5, 1.8], [1]),  # corner, bottomright
        ([(0.5, 0.5), (0.5, 1.4), (1.8, 1.4)], [0.5, 1.4], [2]),  # corner, bottomright
        ([(1.8, 0.5), (1.8, 1.3), (0.9, 1.3)], [0.5, 1.6], [3]),  # corner, topright
        ([(0.5, 1.5), (0.5, 0.2), (1.9, 0.2), (1.9, 1.5)], [0.5, 1.8, 3.5], [0, 1]),
    ],
)
@mock.patch.object(Channels, "get_lines")
def test_embed_channel(get_lines_m, grid2d, channel, lines_s1d, embedded_in, reverse):
    if reverse:
        channel = channel[::-1]
        embedded_in = embedded_in[::-1]
        lines_s1d = pygeos.length(pygeos.linestrings(channel)) - lines_s1d[::-1]
    channels = Channels(
        id=[0],
        calculation_type=EMBEDDED,
        the_geom=[pygeos.linestrings(channel)],
    )

    nodes, lines = embed_channels(
        channels,
        None,
        grid2d.cell_tree,
        count(2),
        count(3),
        connection_node_offset=10,
    )

    assert_array_equal(nodes.embedded_in, embedded_in)
    assert_almost_equal(lines.s1d, lines_s1d)
