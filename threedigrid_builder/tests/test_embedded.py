from numpy.testing import assert_array_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import embed_nodes
from threedigrid_builder.grid import Grid

import pytest


NODE_2D = NodeType.NODE_2D_OPEN_WATER
NODE_1D = NodeType.NODE_1D_STORAGE


@pytest.fixture
def grid2d():
    return Grid(
        nodes=Nodes(
            id=[0, 1],
            node_type=NODE_2D,
            bounds=[(0, 0, 1, 1), (1, 0, 2, 1)],
        ),
        lines=Lines(id=[0], line=[[0, 1]]),
    )


def test_embed_nodes(grid2d):
    EMB = CalculationType.EMBEDDED
    ISO = CalculationType.ISOLATED

    grid = grid2d + Grid(
        nodes=Nodes(
            id=[2, 3, 4, 5],
            node_type=NODE_1D,
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            calculation_type=[EMB, ISO, EMB, EMB],
            coordinates=[(1.0, 0.2), (0.9, 0.9), (1.3, 0.3), (1.6, 0.7)],
            dmax=[1.0, 2.0, 3.0, 4.0],
        ),
        lines=Lines(id=[1, 2, 3], line=[(2, 3), (2, 4), (3, 4)]),
    )
    embedded = embed_nodes(grid)

    assert_array_equal(grid.nodes.id, [0, 1, 2])
    assert_array_equal(grid.nodes.node_type, [NODE_2D] * 2 + [NODE_1D])
    assert_array_equal(grid.lines.line, [(0, 1), (0, 2), (0, 1), (2, 1)])

    assert_array_equal(embedded.id, [0, 1, 2])  # new ids
    assert_array_equal(embedded.calculation_type, EMB)
    assert_array_equal(embedded.embedded_in, [0, 1, 1])


def test_embed_node_outside_2D(grid2d):
    grid = grid2d + Grid(
        nodes=Nodes(
            id=[2, 3],
            node_type=NODE_1D,
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            content_pk=[15, 16],
            calculation_type=CalculationType.EMBEDDED,
            coordinates=[(0.5, 1.2), (0.5, 0.5)],
        ),
        lines=Lines(id=[]),
    )

    with pytest.raises(SchematisationError, match=r".*\[15\] are outside the 2D cells"):
        embed_nodes(grid)


def test_embed_node_two_interconnected(grid2d):
    grid = grid2d + Grid(
        nodes=Nodes(
            id=[2, 3],
            node_type=NODE_1D,
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            content_pk=[4, 5],
            calculation_type=CalculationType.EMBEDDED,
            coordinates=[(0.8, 0.8), (0.5, 0.5)],
        ),
        lines=Lines(id=[2], line=[(2, 3)]),
    )

    with pytest.raises(SchematisationError, match=r".*\[4, 5\] connect to.*"):
        embed_nodes(grid)
