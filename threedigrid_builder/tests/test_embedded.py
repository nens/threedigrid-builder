from itertools import count
from numpy.testing import assert_array_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import Channels
from threedigrid_builder.grid import embed_channel_nodes
from threedigrid_builder.grid import embed_nodes
from threedigrid_builder.grid import Grid, ConnectionNodes

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
    embedded = embed_nodes(grid)

    assert_array_equal(grid.nodes.id, [0, 1, 2, 3, 4])
    assert_array_equal(grid.nodes.node_type, [NODE_2D] * 4 + [NODE_1D])
    assert_array_equal(grid.lines.line, [(0, 1), (0, 4), (0, 1), (4, 1)])

    assert_array_equal(embedded.id, [0, 1, 2])  # new ids
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
        embed_nodes(grid)


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
        embed_nodes(grid)


def test_embed_channel_nodes(grid2d):
    channels = Channels(
        id=[0, 1],
        calculation_type=EMBEDDED,
        the_geom=[
            pygeos.linestrings([(0.1, 0.1), (1.9, 0.1)]),
            pygeos.linestrings([(0.2, 0.2), (1.8, 0.2), (1.8, 1.8)]),
        ],
    )

    nodes = embed_channel_nodes(grid2d.cell_tree, channels, count(2))

    assert_array_equal(nodes.id, [2])
    assert_array_equal(nodes.content_type, ContentType.TYPE_V2_CHANNEL)
    assert_array_equal(nodes.content_pk, [1])
    assert_array_equal(nodes.coordinates, [(1.8, 0.2)])
    assert_array_equal(nodes.calculation_type, EMBEDDED)
    assert_array_equal(nodes.s1d, [1.6])
