from itertools import count

import pytest
from numpy.testing import assert_array_equal

from threedigrid_builder.base import Lines, Nodes
from threedigrid_builder.constants import LineType, NodeType
from threedigrid_builder.grid import Grid
from threedigrid_builder.grid import groundwater as groundwater_module


@pytest.fixture
def grid2d():
    return Grid(
        nodes=Nodes(
            id=[0, 1, 2, 3],
            node_type=NodeType.NODE_2D_OPEN_WATER,
            bounds=[(0, 0, 1, 1), (1, 0, 2, 1), (0, 1, 1, 2), (1, 1, 2, 2)],
            coordinates=[(0.5, 0.5), (1.5, 0.5), (0.5, 1.5), (1.5, 1.5)],
        ),
        lines=Lines(
            id=[0, 1, 2, 3],
            kcu=[
                LineType.LINE_2D_U,
                LineType.LINE_2D_V,
                LineType.LINE_2D_U,
                LineType.LINE_2D_OBSTACLE_V,
            ],
            line=[(0, 1), (0, 2), (1, 3), (2, 3)],
        ),
    )


@pytest.fixture
def grid2d_gw(grid2d):
    grid2d.nodes += groundwater_module.get_nodes(grid2d.nodes, count(4))[0]
    return grid2d


def test_get_nodes(grid2d):
    nodes, id_offset = groundwater_module.get_nodes(grid2d.nodes, count(6))

    assert_array_equal(nodes.id, [6, 7, 8, 9])
    assert_array_equal(nodes.node_type, NodeType.NODE_2D_GROUNDWATER)
    assert_array_equal(nodes.bounds, grid2d.nodes.bounds)
    assert id_offset == 6


def test_get_vertical_lines(grid2d_gw):
    lines = groundwater_module.get_vertical_lines(grid2d_gw.nodes, 4, count(4))

    assert_array_equal(lines.id, [4, 5, 6, 7])
    assert_array_equal(lines.kcu, LineType.LINE_2D_VERTICAL)
    assert_array_equal(lines.line, [(0, 4), (1, 5), (2, 6), (3, 7)])


def test_get_lines(grid2d_gw):
    lines = groundwater_module.get_lines(grid2d_gw.lines, 4, count(4))

    assert_array_equal(lines.id, [4, 5, 6, 7])
    assert_array_equal(lines.kcu, LineType.LINE_2D_GROUNDWATER)
    assert_array_equal(lines.line, [(4, 5), (4, 6), (5, 7), (6, 7)])
