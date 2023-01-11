import itertools
from unittest import mock

import numpy as np
import pygeos
import pytest
from numpy.testing import assert_almost_equal, assert_array_equal
from pygeos.testing import assert_geometries_equal

from threedigrid_builder.base import Lines, Nodes
from threedigrid_builder.constants import CalculationType, ContentType, LineType
from threedigrid_builder.grid import (
    ExchangeLines,
    Lines1D2D,
    Obstacles,
    PotentialBreaches,
)

ISO = CalculationType.ISOLATED
C1 = CalculationType.CONNECTED
C2 = CalculationType.DOUBLE_CONNECTED
CN = ContentType.TYPE_V2_CONNECTION_NODES
CH = ContentType.TYPE_V2_CHANNEL
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
        id=range(1, len(calculation_type) + 1),
        calculation_type=calculation_type,
        content_pk=3,
        content_type=ContentType.TYPE_V2_CHANNEL,
    )
    actual = Lines1D2D.create(nodes, itertools.count())
    assert_array_equal(actual.line[:, 0], -9999)
    assert_array_equal(actual.line[:, 1], expected_line)
    assert_array_equal(actual.content_pk, 3)
    assert_array_equal(actual.content_type, ContentType.TYPE_V2_CHANNEL)


@pytest.mark.parametrize(
    "breach_ids,content_type,expected_content_pk,expected_content_type",
    [
        ((-9999, -9999), CN, 33, CN),
        ((5, -9999), CN, 11, CH),
        ((5, 6), CN, 11, CH),
        ((5, 6), CH, 33, CH),
    ],
)
def test_assign_connection_nodes_to_channels_from_breaches(
    breach_ids, content_type, expected_content_pk, expected_content_type
):
    nodes = Nodes(id=[1], breach_ids=[breach_ids])
    lines = Lines1D2D(
        id=[1], content_type=[content_type], content_pk=[33], line=[(-9999, 1)]
    )
    potential_breaches = PotentialBreaches(id=[5, 6], channel_id=[11, 11])
    lines.assign_connection_nodes_to_channels_from_breaches(nodes, potential_breaches)

    assert_array_equal(lines.content_pk, [expected_content_pk])
    assert_array_equal(lines.content_type, [expected_content_type])


@pytest.mark.parametrize(
    "content_pk,content_type,exc_line_channel_id,expected_content_pk,expected_content_type",
    [
        ([], CH, [], [], []),
        ([11], CH, [], [-9999], [-9999]),
        ([12], CH, [11], [-9999], [-9999]),
        ([12], CH, [12], [1], [EXC]),
        ([12, 12], CH, [12], [1, -9999], [EXC, -9999]),
        ([12, 12], CH, [12, 11], [1, -9999], [EXC, -9999]),
        ([12, 12], CH, [12, 12], [1, 2], [EXC, EXC]),
        ([12], CN, [12], [-9999], [-9999]),
    ],
)
def test_assign_exchange_lines(
    content_pk,
    content_type,
    exc_line_channel_id,
    expected_content_pk,
    expected_content_type,
):
    lines = Lines1D2D(
        id=range(len(content_pk)),
        content_pk=content_pk,
        content_type=content_type,
    )
    exchange_lines = ExchangeLines(
        id=range(1, len(exc_line_channel_id) + 1), channel_id=exc_line_channel_id
    )
    lines.assign_exchange_lines(exchange_lines)

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
    lines.content_type[0] = ContentType.TYPE_V2_BREACH

    lines.assign_kcu(
        mask=np.array([1, 1, 1, 1, 0, 1, 1], dtype=bool),
        is_closed=np.array([0, 0, 1, 1, 0, 1], dtype=bool),
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
@mock.patch.object(Lines1D2D, "assign_ds1d_half_from_obstacles")
def test_assign_dpumax_from_obstacles(
    assign_ds1d_half_from_obstacles,
    assign_dpumax,
):
    obstacles = mock.Mock()
    obstacles.compute_dpumax.return_value = (
        np.array([1.2, np.nan]),
        np.array([1, -9999]),
    )

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

    args, _ = assign_ds1d_half_from_obstacles.call_args

    assert args[0] is obstacles
    assert_array_equal(args[1], [0])
    assert_array_equal(args[2], [1])


@pytest.mark.parametrize(
    "breach_ids,breach_2d_coords,content_pk",
    [
        ([-9999, -9999], [(10, 0)], [-9999]),
        ([1, -9999], [(10, 0)], [1]),
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


def test_assign_breaches_single_connected_err():
    nodes = Nodes(id=[1], breach_ids=[[1, 2]])
    lines = Lines1D2D(id=[1], line=[(-9999, 1)], line_coords=[(10, 0, 0, 0)])
    potential_breaches = PotentialBreaches(
        id=range(1, 3),
    )
    with pytest.raises(ValueError):
        lines.assign_breaches(nodes, potential_breaches)


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


@pytest.mark.parametrize(
    "obstacle_geom,line_coords,expected",
    [
        ([[5, 0], [5, 10]], [0, 5, 10, 5], 5.0),
        ([[5, 0], [5, 10]], [2, 5, 10, 5], 3.0),
        ([[5, 0], [5, 10]], [0, 0, 10, 0], 5.0),
        ([[5, 1], [5, 10]], [0, 0, 10, 0], 5.0),
        ([[5, 3], [5, 10]], [5, 0, 5, 10], 3.0),
    ],
)
def test_assign_ds1d_half_from_obstacles(obstacle_geom, line_coords, expected):
    obstacles = Obstacles(id=[1], the_geom=[pygeos.linestrings(obstacle_geom)])
    lines = Lines1D2D(id=range(1), line_coords=[line_coords])
    lines.fix_line_geometries()

    lines.assign_ds1d_half_from_obstacles(obstacles, np.array([0]), np.array(([0])))

    assert_almost_equal(lines.ds1d_half, [expected])


def test_assign_ds1d():
    nodes = Nodes(id=[1, 2, 3], bounds=[(0, 0, 8, 8), (8, 0, 10, 2), (8, 2, 10, 4)])
    lines = Lines1D2D(
        id=range(4),
        line=[[1, -9999], [1, -9999], [2, -9999], [3, -9999]],
        ds1d=[np.nan, np.nan, np.nan, 5.2],
    )
    lines.assign_ds1d(nodes)

    assert_array_equal(lines.ds1d, [8.0, 8.0, 2.0, 5.2])


def test_assign_ds1d_half():
    lines = Lines1D2D(
        id=range(2),
        ds1d_half=[np.nan, 5.0],
        line_geometries=[pygeos.linestrings([[0, 0], [0, 9]]), None],
    )
    lines.assign_ds1d_half()

    assert_array_equal(lines.ds1d_half, [4.5, 5.0])


@pytest.fixture
def threeway_junction():
    nodes = Nodes(
        id=range(1, 5),
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        content_pk=[1, 2, 3, 4],
    )
    lines = Lines(
        id=range(3),
        line=[(1, 2), (1, 3), (4, 1)],
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=[11, 12, 13],
    )
    return nodes, lines


@pytest.fixture
def potential_breaches():
    return PotentialBreaches(
        id=[5, 6, 7],
        channel_id=[11, 12, 13],
    )


def test_channel_mapping():
    exchange_lines = ExchangeLines(id=range(1, 6), channel_id=[15, 2, 3, 15, 3])
    actual = exchange_lines.channel_mapping

    assert_array_equal(actual.id, [2, 3, 15])
    assert_array_equal(actual.exchange_line_id, [2, 3, 1])
    assert_array_equal(actual.secondary_exchange_line_id, [-9999, 5, 4])


@pytest.mark.parametrize(
    "channel_id,expected,is_primary",
    [
        (np.array([]), [], True),
        (np.array([1]), [-9999], True),
        (np.array([2]), [2], True),
        (np.array([2]), [-9999], False),
        (np.array([3]), [3], True),
        (np.array([3]), [5], False),
    ],
)
def test_get_for_channel_id(channel_id, expected, is_primary):
    exchange_lines = ExchangeLines(id=range(1, 6), channel_id=[15, 2, 3, 15, 3])
    actual = exchange_lines.get_for_channel_id(channel_id, is_primary)

    assert_array_equal(actual, expected)


@pytest.mark.parametrize(
    "kcu,expected",
    [
        ([ISO, ISO, ISO], 33),
        ([C1, ISO, ISO], 11),
        ([ISO, ISO, C1], 13),
        ([ISO, C1, C1], 12),
        ([C2, ISO, ISO], 11),
        ([ISO, ISO, C2], 13),
        ([ISO, C2, C2], 12),
        ([C1, C2, C2], 12),
    ],
)
@pytest.mark.parametrize("is_double_connected", [True, False])
def test_assign_connection_nodes_to_channels_from_lines(
    threeway_junction, kcu, expected, is_double_connected
):
    threeway_junction[1].kcu[:] = kcu
    n = 2 if is_double_connected else 1
    lines_1d2d = Lines1D2D(
        id=range(n), content_type=CN, content_pk=33, line=[(-9999, 1)] * n
    )
    lines_1d2d.assign_connection_nodes_to_channels_from_lines(*threeway_junction)

    assert_array_equal(lines_1d2d.content_pk, expected)
    assert_array_equal(lines_1d2d.content_type, CN if expected == 33 else CH)


@mock.patch.object(Lines1D2D, "get_velocity_points")
def test_output_breaches(get_velocity_points):
    lines = Lines1D2D(
        id=range(1, 5),
        content_type=[-9999, BREACH, EXC, BREACH],
        content_pk=[-9999, 1, 2, 3],
    )
    potential_breaches = PotentialBreaches(
        id=[1, 2, 3],
        code=["a", "b", "c"],
        display_name=["aa", "bb", "cc"],
        levee_material=[1, 2, -9999],
        maximum_breach_depth=[np.nan, 1.3, 1.4],
    )

    get_velocity_points.return_value = pygeos.points([[0, 0], [1, 1], [2, 2]])
    actual = lines.output_breaches(potential_breaches)

    assert_array_equal(actual.id, [0, 1, 2])
    assert_array_equal(actual.line_id, [2, 4, 3])
    assert_array_equal(actual.content_pk, [1, 3, -9999])

    assert_array_equal(actual.maximum_breach_depth, [np.nan, 1.4, np.nan])
    assert_array_equal(actual.levee_material, [1, -9999, -9999])

    assert_geometries_equal(actual.the_geom, get_velocity_points.return_value)

    assert_array_equal(actual.code, ["a", "c", None])
    assert_array_equal(actual.display_name, ["aa", "cc", None])

    get_velocity_points.assert_called_once()
