import itertools
from unittest import mock

import numpy as np
import pygeos
import pytest
from numpy.testing import assert_array_equal

from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType, ContentType, LineType
from threedigrid_builder.grid import ExchangeLines, Lines1D2D, PotentialBreaches

ISO = CalculationType.ISOLATED
C1 = CalculationType.CONNECTED
C2 = CalculationType.DOUBLE_CONNECTED
EXC = ContentType.TYPE_V2_EXCHANGE_LINE
BREACH = ContentType.TYPE_V2_BREACH


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


def test_assign_2d_side():
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
    lines.assign_2d_side(nodes, exchange_lines)

    assert_array_equal(lines.line_coords[:, :2], [[2, 0], [2, 10], [4, 0], [6, 5]])


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
        line_coords=[x + (np.nan, np.nan) for x in side_2d_coordinates],
    )

    lines.assign_2d_node(cell_tree)

    assert_array_equal(lines.line[:, 0], expected_2d_node_id)
    assert_array_equal(lines.line_coords, np.nan)  # is cleared


def test_assign_kcu():
    lines = Lines1D2D(id=range(7), line=[[-9999] + [x] for x in [1, 1, 2, 2, 3, 4, 5]])

    lines.assign_kcu(
        np.array([1, 1, 1, 1, 0, 1, 1], dtype=bool),
        np.array([0, 0, 1, 1, 0, 1], dtype=bool),
    )

    assert_array_equal(
        lines.kcu,
        [
            LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_CLOSED,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_CLOSED,
            -9999,
            LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER,
            LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED,
        ],
    )


@pytest.mark.parametrize(
    "existing,mask,expected",
    [
        ([np.nan, np.nan, np.nan], [1, 0, 1], [3.0, np.nan, 6.0]),
        ([5.0, 5.0, 5.0], [1, 0, 1], [5.0, 5.0, 5.0]),
        ([5.0, np.nan, np.nan], [1, 0, 1], [5.0, np.nan, 6.0]),
    ],
)
def test_assign_dpumax(existing, mask, expected):
    lines = Lines1D2D(id=range(3), dpumax=existing)

    lines.assign_dpumax(
        np.array(mask, dtype=bool),
        np.array([3.0, 6.0]),
    )

    assert_array_equal(lines.dpumax, expected)


def test_get_1d_node_idx():
    nodes = Nodes(
        id=[1, 2, 3],
        coordinates=[(2, 5), (4, 5), (6, 5)],
    )
    lines = Lines1D2D(
        id=range(4),
        line=[[-9999, 1], [-9999, 1], [-9999, 2], [-9999, 3]],
    )
    actual = lines.get_1d_node_idx(nodes)

    assert_array_equal(actual, [0, 0, 1, 2])


@mock.patch.object(Lines1D2D, "assign_dpumax")
def test_assign_dpumax_from_exchange_lines(assign_dpumax):
    lines = Lines1D2D(id=range(3), content_pk=[1, 2, 3], content_type=[EXC, -9999, EXC])
    exchange_lines = ExchangeLines(id=[1, 2, 3], exchange_level=[1.2, 2.3, np.nan])

    lines.assign_dpumax_from_exchange_lines(exchange_lines)

    (actual_mask, actual_dpumax), _ = assign_dpumax.call_args

    assert_array_equal(actual_mask, [1, 0, 1])
    assert_array_equal(actual_dpumax, [1.2, np.nan])


@mock.patch.object(Lines1D2D, "assign_dpumax")
def test_assign_dpumax_from_breaches(assign_dpumax):
    lines = Lines1D2D(
        id=range(3), content_pk=[1, 2, 3], content_type=[BREACH, -9999, BREACH]
    )
    breaches = PotentialBreaches(id=[1, 2, 3], exchange_level=[1.2, 2.3, np.nan])

    lines.assign_dpumax_from_breaches(breaches)

    (actual_mask, actual_dpumax), _ = assign_dpumax.call_args

    assert_array_equal(actual_mask, [1, 0, 1])
    assert_array_equal(actual_dpumax, [1.2, np.nan])


@mock.patch.object(Lines1D2D, "assign_dpumax")
def test_assign_dpumax_from_obstacles(assign_dpumax):
    obstacles = mock.Mock()
    obstacles.compute_dpumax.return_value = (np.array([1.2, np.nan]), None)

    lines = Lines1D2D(id=range(3), content_type=[EXC, -9999, EXC])

    lines.assign_dpumax_from_obstacles(obstacles)

    args, kwargs = obstacles.compute_dpumax.call_args
    assert len(args) == 1
    assert args[0] is lines
    assert len(kwargs) == 1
    assert_array_equal(kwargs["where"], [0, 2])

    (actual_mask, actual_dpumax), _ = assign_dpumax.call_args

    assert_array_equal(actual_mask, [1, 0, 1])
    assert_array_equal(actual_dpumax, [1.2, np.nan])


@pytest.mark.parametrize(
    "breach_ids,breach_2d_coords,content_pk",
    [
        ([-9999, -9999], [(10, 0)], [-9999]),
        ([1, -9999], [(10, 0)], [1]),
        ([1, 2], [(10, 0), (10, 1)], [1]),
        ([1, 2], [(10, 1), (10, 0)], [2]),
        ([1, 2], [(10, 1), (10, -1)], [1]),
    ],
)
def test_assign_breaches_single_connected(breach_ids, breach_2d_coords, content_pk):
    nodes = Nodes(id=[1], breach_ids=[breach_ids])
    lines = Lines1D2D(id=[1], line=[(-9999, 1)], line_coords=[(10, 0, 0, 0)])
    potential_breaches = PotentialBreaches(
        id=range(1, len(breach_2d_coords) + 1),
        the_geom=pygeos.linestrings([[(0, 0), x] for x in breach_2d_coords]),
    )
    lines.assign_breaches(nodes, potential_breaches)

    assert_array_equal(
        lines.content_type, np.where(np.array(content_pk) != -9999, BREACH, -9999)
    )
    assert_array_equal(lines.content_pk, content_pk)


@pytest.mark.parametrize(
    "breach_ids,breach_2d_coords,content_pk",
    [
        ([-9999, -9999], [(10, 0)], [-9999, -9999]),
        ([1, -9999], [(10, 0)], [1, -9999]),
        ([1, -9999], [(10, 10)], [-9999, 1]),
        ([1, -9999], [(10, 5)], [1, -9999]),
        ([1, 2], [(10, 0), (10, 10)], [1, 2]),
        ([1, 2], [(10, 10), (10, 0)], [2, 1]),
        ([1, 2], [(10, 1), (10, 9)], [1, 2]),
        ([1, 2], [(10, 9), (10, 1)], [2, 1]),
        ([1, 2], [(10, 5), (10, 5)], [1, 2]),
    ],
)
def test_assign_breaches_double_connected(breach_ids, breach_2d_coords, content_pk):
    nodes = Nodes(id=[1], breach_ids=[breach_ids])
    lines = Lines1D2D(
        id=[1, 2],
        line=[(-9999, 1), (-9999, 1)],
        line_coords=[(10, 0, 0, 0), (10, 10, 0, 0)],
    )
    potential_breaches = PotentialBreaches(
        id=range(1, len(breach_2d_coords) + 1),
        the_geom=pygeos.linestrings([[(0, 0), x] for x in breach_2d_coords]),
    )
    lines.assign_breaches(nodes, potential_breaches)

    assert_array_equal(
        lines.content_type, np.where(np.array(content_pk) != -9999, BREACH, -9999)
    )
    assert_array_equal(lines.content_pk, content_pk)


def test_assign_breaches_multiple():
    nodes = Nodes(
        id=[1, 2, 3],
        coordinates=[(2, 5), (4, 5), (6, 5)],
        breach_ids=[(1, 2), (-9999, -9999), (4, -9999)],
    )
    lines = Lines1D2D(
        id=range(4),
        line=[[-9999, 1], [-9999, 1], [-9999, 2], [-9999, 3]],
        line_coords=[x + [0, 0] for x in [[2, 0], [2, 10], [4, 0], [6, 5]]],
        content_type=[EXC, EXC, EXC, -9999],
        content_pk=[0, 1, 2, 3],
    )
    potential_breaches = PotentialBreaches(
        id=[1, 2, 3, 4],
        the_geom=pygeos.linestrings(
            [[(0, 0), x] for x in [(2, 9), (2, 1), (6, 10), (6, 1)]]
        ),
        code=["a", "b", "c", "d"],
    )
    lines.assign_breaches(nodes, potential_breaches)

    # also check the other fields that are set by assign_breaches
    assert_array_equal(lines.line_coords[:, :2], [(2, 1), (2, 9), (4, 0), (6, 1)])
    assert_array_equal(lines.content_type, [BREACH, BREACH, EXC, BREACH])
    assert_array_equal(lines.content_pk, [2, 1, 2, 4])
