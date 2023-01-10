import logging

import numpy as np
import pygeos
import pytest
from numpy.testing import assert_almost_equal, assert_equal

from threedigrid_builder.base import Lines, Nodes
from threedigrid_builder.constants import CalculationType, ContentType, Material
from threedigrid_builder.grid import (
    Channels,
    ConnectedPoints,
    Levees,
    Obstacles,
    PotentialBreaches,
    PotentialBreachPoints,
)

CH = ContentType.TYPE_V2_CHANNEL
CP = ContentType.TYPE_V2_ADDED_CALCULATION_POINT
BREACH = ContentType.TYPE_V2_BREACH


@pytest.fixture
def levees():
    return Levees(
        id=[1, 2, 3],
        the_geom=[
            pygeos.linestrings([[0, 0], [10, 0], [10, 10]]),
            pygeos.linestrings([[5, 5], [8, 5]]),
            pygeos.linestrings([[0, 0], [1, 1]]),
        ],
        crest_level=[2.0, 4.0, np.nan],
        max_breach_depth=[4.0, 2.0, np.nan],
        material=[Material.CLAY, Material.SAND, -9999],
    )


@pytest.fixture
def connected_points():
    return ConnectedPoints(
        id=[0, 1, 2, 3],
        levee_id=[1, 1, 2, 3],
    )


@pytest.fixture
def lines():
    return Lines(
        id=[0, 1, 2, 3, 4, 5],
        line_geometries=[
            None,
            None,
            pygeos.linestrings([[5, 5], [15, 5]]),  # crosses levee 1 at [10, 5]
            pygeos.linestrings([[3, 2], [4, 1], [5, 2]]),  # (closest on levee: [4, 0])
            pygeos.linestrings([[6, 5], [10, 5]]),  # tangent to levee 2
            None,
        ],
        content_type=[CH, -9999, CP, CP, CP, CP],
        content_pk=[1, -9999, 0, 1, 2, 3],
    )


def test_get_breaches(connected_points, lines, levees):
    breaches = connected_points.get_breaches(lines, levees)

    assert isinstance(breaches, PotentialBreaches)
    assert len(breaches) == 3

    assert_equal(lines.content_type, [CH, -9999, BREACH, BREACH, BREACH, CP])
    assert_equal(lines.content_pk, [1, -9999, 0, 1, 2, 3])
    assert_almost_equal(
        lines.ds1d_half, [np.nan, np.nan, 5.0, 1.4142, 0.0, np.nan], decimal=4
    )

    assert_equal(breaches.id, [0, 1, 2])
    assert_almost_equal(breaches.maximum_breach_depth, [4, 4, 2])
    assert_equal(breaches.levee_material, [Material.CLAY, Material.CLAY, Material.SAND])


def test_no_breaches(connected_points, lines, levees):
    connected_points.levee_id[:] = -9999
    breaches = connected_points.get_breaches(lines, levees)

    assert isinstance(breaches, PotentialBreaches)
    assert len(breaches) == 0


def test_potential_breach_sides():
    potential_breaches = PotentialBreaches(
        id=[1],
        channel_id=[1],
        the_geom=[pygeos.linestrings([[0, 0], [0, 1]])],
    )
    assert_almost_equal(pygeos.get_coordinates(potential_breaches.side_1d), [[0, 0]])
    assert_almost_equal(pygeos.get_coordinates(potential_breaches.side_2d), [[0, 1]])


def test_potential_breach_merge(caplog):
    breaches = PotentialBreachPoints(
        linestrings=Channels(
            id=[0, 1],
            the_geom=pygeos.linestrings([[[0, 0], [10, 0]], [[0, 0], [0, 10]]]),
        ).linestrings,
        id=range(7),
        s1d=[0.00001, 2.0, 5.0, 5.00001, 9.99999, 4, 4.00001],
        linestring_idx=[0, 0, 0, 0, 0, 1, 1],
        content_pk=range(1, 8),
    )

    with caplog.at_level(logging.WARNING):
        actual = breaches.merge()

    assert len(caplog.messages) == 0

    assert_almost_equal(actual.s1d, [0.0, 2.0, 5.0, 10.0, 4.0])
    assert_almost_equal(actual.linestring_idx, [0, 0, 0, 0, 1])
    assert_almost_equal(actual.content_pk, [1, 2, 3, 5, 6])
    assert_almost_equal(actual.secondary_content_pk, [-9999, -9999, 4, -9999, 7])


def test_potential_breach_merge_empty():
    actual = PotentialBreachPoints.empty(
        linestrings=Channels(
            id=[0, 1],
            the_geom=pygeos.linestrings([[[0, 0], [10, 0]], [[0, 0], [0, 10]]]),
        ).linestrings
    ).merge()

    assert len(actual) == 0


def test_levees_merge_into_obstacles(levees):
    obstacles = Obstacles(id=[2, 3], crest_level=[5.0, np.nan])
    actual = levees.merge_into_obstacles(obstacles)

    assert isinstance(actual, Obstacles)
    assert actual is not obstacles

    # levees are added, with other ids
    assert_equal(actual.id, [2, 3, 4, 5, 6])
    assert_equal(actual.crest_level, [5.0, np.nan, 2.0, 4.0, np.nan])


def test_levees_merge_into_obstacles_no_obstacles(levees):
    obstacles = Obstacles(id=[])
    actual = levees.merge_into_obstacles(obstacles)

    assert isinstance(actual, Obstacles)
    assert actual is not obstacles

    assert_equal(actual.id, [1, 2, 3])
    assert_equal(actual.crest_level, [2.0, 4.0, np.nan])


@pytest.fixture
def threeway_junction():
    nodes = Nodes(
        id=range(1, 5),
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        content_pk=[1, 2, 3, 4],
    )
    linestrings = Channels(
        id=[11, 12, 13],
        the_geom=pygeos.linestrings([[[0, 0], [0, 1]]] * 3),
    ).linestrings
    lines = Lines(
        id=range(3),
        line=[(1, 2), (1, 3), (4, 1)],
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=[11, 12, 13],
    )
    return nodes, lines, linestrings


@pytest.mark.parametrize(
    "breach_points,expected",
    [
        (  # no breach points: nothing assigned
            {"id_1": [], "id_2": [], "s1d": [], "ch_idx": []},
            (-9999, -9999),
        ),
        (  # breach point at start, single connected
            {"id_1": [5], "id_2": [-9999], "s1d": [0.0], "ch_idx": [0]},
            (5, -9999),
        ),
        (  # breach point at end, single connected (linestring 2 is swapped)
            {"id_1": [7], "id_2": [-9999], "s1d": [1.0], "ch_idx": [2]},
            (7, -9999),
        ),
        (  # breach point at start, far side
            {"id_1": [5], "id_2": [-9999], "s1d": [1.0], "ch_idx": [0]},
            (-9999, -9999),
        ),
        (  # breach point at end, far side (linestring 2 is swapped)
            {"id_1": [7], "id_2": [-9999], "s1d": [0.0], "ch_idx": [2]},
            (-9999, -9999),
        ),
        (  # breach point at start, double connected
            {"id_1": [5], "id_2": [6], "s1d": [0.0], "ch_idx": [0]},
            (5, 6),
        ),
        (  # breach point at end, double connected (linestring 2 is swapped)
            {"id_1": [7], "id_2": [8], "s1d": [1.0], "ch_idx": [2]},
            (7, 8),
        ),
        (  # two breach points, single connected
            {"id_1": [5, 6], "id_2": -9999, "s1d": [0.0, 0.0], "ch_idx": [0, 1]},
            (5, -9999),
        ),
        (  # two breach points, single connected (orders by ch_idx)
            {"id_1": [5, 6], "id_2": -9999, "s1d": [0.0, 0.0], "ch_idx": [1, 0]},
            (6, -9999),
        ),
        (  # two breach points, one double connected (gets priority)
            {"id_1": [5, 6], "id_2": [-9999, 8], "s1d": [0.0, 0.0], "ch_idx": [0, 1]},
            (6, 8),
        ),
    ],
)
def test_assign_to_connection_nodes(threeway_junction, breach_points, expected):
    """Set the breach ids on a connection node with 3 channels"""
    nodes, lines, linestrings = threeway_junction
    breach_points = PotentialBreachPoints(
        linestrings=linestrings,
        id=range(len(breach_points["id_1"])),
        content_pk=breach_points["id_1"],
        secondary_content_pk=breach_points["id_2"],
        s1d=breach_points["s1d"],
        linestring_idx=breach_points["ch_idx"],
    )
    breach_points.assign_to_connection_nodes(nodes, lines)
    assert_equal(nodes.breach_ids[0], expected)


@pytest.mark.parametrize(
    "calculation_type,breach_ids,expected,messages",
    [
        (CalculationType.ISOLATED, [-9999, -9999], [-9999, -9999], []),
        (CalculationType.CONNECTED, [1, -9999], [1, -9999], []),
        (CalculationType.DOUBLE_CONNECTED, [1, 2], [1, 2], []),
        (
            CalculationType.ISOLATED,
            [1, -9999],
            [-9999, -9999],
            [
                "The following objects have potential breaches, but are not (double) connected: channels [11].",
                "The following potential breaches will be ignored: [1].",
            ],
        ),
        (
            CalculationType.ISOLATED,
            [1, 2],
            [-9999, -9999],
            [
                "The following objects have potential breaches, but are not (double) connected: channels [11].",
                "The following potential breaches will be ignored: [1, 2].",
            ],
        ),
        (
            CalculationType.CONNECTED,
            [1, 2],
            [1, -9999],
            [
                "The following objects have two potential breaches at the same position, but are not double connected: channels [11].",
                "The following potential breaches will be ignored: [2].",
            ],
        ),
    ],
)
def test_match_breach_ids_with_calculation_types(
    caplog, calculation_type, breach_ids, expected, messages
):
    nodes = Nodes(
        id=[1],
        calculation_type=calculation_type,
        breach_ids=[breach_ids],
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=[11],
    )
    nodes.breach_ids[0, :] = breach_ids

    with caplog.at_level(logging.WARNING):
        PotentialBreachPoints.match_breach_ids_with_calculation_types(nodes)

    assert caplog.messages == messages
    assert_equal(nodes.breach_ids[0], expected)
