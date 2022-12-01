import itertools

import numpy as np
import pygeos
import pytest
from numpy.testing import assert_array_equal

from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType, ContentType
from threedigrid_builder.grid import ExchangeLines, Lines1D2D

ISO = CalculationType.ISOLATED
C1 = CalculationType.CONNECTED
C2 = CalculationType.DOUBLE_CONNECTED
EXC = ContentType.TYPE_V2_EXCHANGE_LINE


@pytest.mark.parametrize(
    "calculation_type,expected_line",
    [
        ([], []),
        ([ISO], []),
        ([C1], [1]),
        ([C2], [1, 1]),
        ([C1, ISO, C2], [1, 3, 3]),
        ([C2, ISO, C1], [1, 1, 3]),
    ],
)
def test_create_lines_1d2d(calculation_type, expected_line):
    nodes = Nodes(
        id=range(1, len(calculation_type) + 1), calculation_type=calculation_type
    )
    actual = Lines1D2D.create(nodes, itertools.count())
    assert_array_equal(actual.line[:, 0], -9999)
    assert_array_equal(actual.line[:, 1], expected_line)


@pytest.mark.parametrize(
    "line_node_id,exc_line_channel_id,expected_content_pk,expected_content_type",
    [
        ([], [], [], []),
        ([2], [], [-9999], [-9999]),
        ([2], [11], [-9999], [-9999]),
        ([2], [12], [1], [EXC]),
        ([2, 2], [12], [1, -9999], [EXC, -9999]),
        ([2, 2], [12, 11], [1, -9999], [EXC, -9999]),
        ([2, 2], [12, 12], [1, 2], [EXC, EXC]),
    ],
)
def test_assign_exchange_lines(
    line_node_id, exc_line_channel_id, expected_content_pk, expected_content_type
):
    nodes = Nodes(
        id=[2],
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=[12],
    )
    line = np.full((len(line_node_id), 2), -9999, dtype=np.int32)
    line[:, 1] = line_node_id
    lines = Lines1D2D(id=range(len(line_node_id)), line=line)
    exchange_lines = ExchangeLines(
        id=range(1, len(exc_line_channel_id) + 1), channel_id=exc_line_channel_id
    )

    lines.assign_exchange_lines(nodes, exchange_lines)

    assert_array_equal(lines.content_pk, expected_content_pk)
    assert_array_equal(lines.content_type, expected_content_type)


def test_assign_line_coords():
    nodes = Nodes(
        id=[1, 2, 3],
        coordinates=[(2, 5), (4, 5), (6, 5)],
    )
    exchange_lines = ExchangeLines(
        id=[1, 2],
        the_geom=pygeos.linestrings([[[0, 0], [10, 0]], [[0, 10], [10, 10]]]),
    )
    lines = Lines1D2D(
        id=range(4),
        line=[[-9999, 1], [-9999, 1], [-9999, 2], [-9999, 3]],
        content_pk=[1, 2, 1, -9999],
        content_type=[EXC, EXC, EXC, -9999],
    )
    lines.assign_line_coords(nodes, exchange_lines)

    assert_array_equal(
        lines.line_coords, [[2, 0, 2, 5], [2, 10, 2, 5], [4, 0, 4, 5], [6, 5, 6, 5]]
    )


@pytest.fixture
def cell_tree():
    return pygeos.STRtree([pygeos.box(0, 0, 1, 1), pygeos.box(1, 0, 3, 2)])


@pytest.mark.parametrize(
    "side_2d_coordinates,expected_2d_node_id",
    [
        ([(0.5, 0.5)], [0]),  # first cell, center
        ([(0, 0.5)], [0]),  # first cell, left edge
        ([(0.5, 1)], [0]),  # first cell, top edge
        ([(0.5, 0)], [0]),  # first cell, bottom edge
        ([(0, 1)], [0]),  # first cell, topleft corner
        ([(0, 0)], [0]),  # first cell, bottomleft corner
        ([(2, 1)], [1]),  # second cell, center
        ([(3, 1)], [1]),  # second cell, right edge
        ([(2, 2)], [1]),  # second cell, top edge
        ([(2, 0)], [1]),  # second cell, bottom edge
        ([(3, 2)], [1]),  # second cell, topright corner
        ([(3, 0)], [1]),  # second cell, bottomright corner
        ([(1, 1)], [0]),  # edge between: top corner
        ([(1, 0)], [0]),  # edge between: bottom corner
        ([(1, 0.5)], [0]),  # edge between: middle
        ([(0.5, 0.5), (0.5, 0.9)], [0, 0]),  # two cells, same
        ([(0.5, 0.5), (2, 1)], [0, 1]),  # two cells, different
        ([(3, 3)], [-9999]),  # out of bounds
    ],
)
def test_assign_2d_node(side_2d_coordinates, expected_2d_node_id, cell_tree):
    lines = Lines1D2D(
        id=range(len(side_2d_coordinates)),
        line_coords=[x + (0, 0) for x in side_2d_coordinates],
    )

    lines.assign_2d_node(cell_tree)

    assert_array_equal(lines.line[:, 0], expected_2d_node_id)
