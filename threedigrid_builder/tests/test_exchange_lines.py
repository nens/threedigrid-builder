from unittest import mock

import numpy as np
import pygeos
import pytest
from numpy.testing import assert_array_equal
from pygeos.testing import assert_geometries_equal

from threedigrid_builder.base import Lines, Nodes
from threedigrid_builder.constants import CalculationType, ContentType, NodeType
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import ExchangeLines, Grid


@pytest.mark.parametrize(
    "content_type",
    [
        ContentType.TYPE_V2_CONNECTION_NODES,
        ContentType.TYPE_V2_PIPE,
        ContentType.TYPE_V2_CULVERT,
    ],
)
def test_get_line_mappings_non_channels(content_type):
    nodes = Nodes(
        id=[1],
        content_type=content_type,
        content_pk=[2],
        calculation_type=[CalculationType.CONNECTED],
    )
    exchange_lines = ExchangeLines(id=[1], channel_id=[2])
    line_node_idx, line_exc_id = exchange_lines.get_line_mappings(nodes)

    assert_array_equal(line_node_idx, [0])
    assert_array_equal(line_exc_id, [-9999])


ISO = CalculationType.ISOLATED
C1 = CalculationType.CONNECTED
C2 = CalculationType.DOUBLE_CONNECTED


@pytest.mark.parametrize(
    "node_calc_type,node_content_pk,channel_id,expected_node_idx,expected_exc_id",
    [
        ([], [], [], [], []),
        ([ISO], [1], [], [], []),
        ([C1], [1], [], [0], [-9999]),
        ([C1], [1], [1], [0], [1]),
        ([C1], [2], [1], [0], [-9999]),
        ([C2], [1], [], [0, 0], [-9999, -9999]),
        ([C2], [1], [1], [0, 0], [1, -9999]),
        ([C2], [2], [1], [0, 0], [-9999, -9999]),
        ([C2], [1], [1, 1], [0, 0], [1, 2]),
        ([C2], [1], [1, 2], [0, 0], [1, -9999]),
        ([C1, ISO, C2], [1, 2, 3], [1, 3, 3], [0, 2, 2], [1, 2, 3]),
    ],
)
def test_get_line_mappings_channels(
    node_calc_type, node_content_pk, channel_id, expected_node_idx, expected_exc_id
):
    nodes = Nodes(
        id=range(len(node_content_pk)),
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=node_content_pk,
        calculation_type=node_calc_type,
    )
    exchange_lines = ExchangeLines(
        id=range(1, len(channel_id) + 1), channel_id=channel_id
    )

    line_node_idx, line_exc_id = exchange_lines.get_line_mappings(nodes)

    assert_array_equal(line_node_idx, expected_node_idx)
    assert_array_equal(line_exc_id, expected_exc_id)


def test_get_line_mappings_too_many_exchange_lines():
    exchange_lines = ExchangeLines(id=[1, 2, 3], channel_id=[1, 1, 1])

    with pytest.raises(
        SchematisationError, match=r"Channels \[1\] have too many exchange lines"
    ):
        exchange_lines.get_line_mappings(Nodes(id=[]))


def test_get_2d_side_points():
    nodes = Nodes(
        id=[0, 1, 2],
        coordinates=[(2, 5), (4, 5), (6, 5)],
    )
    exchange_lines = ExchangeLines(
        id=[1, 2],
        the_geom=pygeos.linestrings([[[0, 0], [10, 0]], [[0, 10], [10, 10]]]),
    )
    # nodes have respectively: 2 exchange lines, 1 exchange line, no exchange lines
    node_idx = np.array([0, 0, 1, 2])
    exc_ids = np.array([1, 2, 1, -9999])

    expected_geometries = pygeos.points([[2, 0], [2, 10], [4, 0], [6, 5]])

    with mock.patch.object(
        exchange_lines, "get_line_mappings", return_value=(node_idx, exc_ids)
    ):
        actual = exchange_lines.get_2d_side_points(nodes)

        assert actual[0] is node_idx
        assert actual[1] is exc_ids
        actual_geometries = actual[2]

    assert_geometries_equal(actual_geometries, expected_geometries)


@pytest.fixture
def grid2d():
    return Grid(
        nodes=Nodes(
            id=[0, 1],
            node_type=NodeType.NODE_2D_OPEN_WATER,
            bounds=[(0, 0, 1, 1), (1, 0, 3, 2)],
        ),
        lines=Lines(id=[0]),
    )


@pytest.mark.parametrize(
    "side_2d_coordinates,expected_lines",
    [
        ([(0.5, 0.5)], [(0, 2)]),  # first cell, center
        ([(0, 0.5)], [(0, 2)]),  # first cell, left edge
        ([(0.5, 1)], [(0, 2)]),  # first cell, top edge
        ([(0.5, 0)], [(0, 2)]),  # first cell, bottom edge
        ([(0, 1)], [(0, 2)]),  # first cell, topleft corner
        ([(0, 0)], [(0, 2)]),  # first cell, bottomleft corner
        ([(2, 1)], [(1, 2)]),  # second cell, center
        ([(3, 1)], [(1, 2)]),  # second cell, right edge
        ([(2, 2)], [(1, 2)]),  # second cell, top edge
        ([(2, 0)], [(1, 2)]),  # second cell, bottom edge
        ([(3, 2)], [(1, 2)]),  # second cell, topright corner
        ([(3, 0)], [(1, 2)]),  # second cell, bottomright corner
        ([(1, 1)], [(0, 2)]),  # edge between: top corner
        ([(1, 0)], [(0, 2)]),  # edge between: bottom corner
        ([(1, 0.5)], [(0, 2)]),  # edge between: middle
        ([(0.5, 0.5), (0.5, 0.9)], [(0, 2), (0, 3)]),  # two cells, same
        ([(0.5, 0.5), (2, 1)], [(0, 2), (1, 3)]),  # two cells, different
    ],
)
def test_get_lines_node_idx(side_2d_coordinates, expected_lines, grid2d):
    grid2d.nodes += Nodes(
        id=[7, 8][: len(side_2d_coordinates)],
        coordinates=side_2d_coordinates,
    )
    exchange_lines = ExchangeLines(id=[])
    node_idx = np.array([2, 3][: len(side_2d_coordinates)])
    exc_ids = np.array([-9999, -9999][: len(side_2d_coordinates)])
    side_2d_geometries = pygeos.points(side_2d_coordinates)

    with mock.patch.object(
        exchange_lines,
        "get_2d_side_points",
        return_value=(node_idx, exc_ids, side_2d_geometries),
    ):
        idx_1d, idx_2d, exc_ids_ = exchange_lines.get_lines_node_idx(
            grid2d.nodes, grid2d.cell_tree
        )

        assert exc_ids is exc_ids_

    expected_idx_1d, expected_idx_2d = np.array(expected_lines).T
    assert_array_equal(idx_1d, expected_idx_1d)
    assert_array_equal(idx_2d, expected_idx_2d)
