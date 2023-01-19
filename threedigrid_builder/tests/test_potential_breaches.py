import logging

import pytest
import shapely
from numpy.testing import assert_almost_equal, assert_equal

from threedigrid_builder.base import Lines, Nodes
from threedigrid_builder.constants import CalculationType, ContentType
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import Channels, PotentialBreaches, PotentialBreachPoints


def test_potential_breach_sides():
    potential_breaches = PotentialBreaches(
        id=[1],
        channel_id=[1],
        the_geom=[shapely.linestrings([[0, 0], [0, 1]])],
    )
    assert_almost_equal(shapely.get_coordinates(potential_breaches.side_1d), [[0, 0]])
    assert_almost_equal(shapely.get_coordinates(potential_breaches.side_2d), [[0, 1]])


def test_potential_breach_merge(caplog):
    breaches = PotentialBreachPoints(
        linestrings=Channels(
            id=[0, 1],
            the_geom=shapely.linestrings([[[0, 0], [10, 0]], [[0, 0], [0, 10]]]),
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
            the_geom=shapely.linestrings([[[0, 0], [10, 0]], [[0, 0], [0, 10]]]),
        ).linestrings
    ).merge()

    assert len(actual) == 0


@pytest.fixture
def threeway_junction():
    nodes = Nodes(
        id=range(1, 5),
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        content_pk=[1, 2, 3, 4],
    )
    linestrings = Channels(
        id=[11, 12, 13],
        the_geom=shapely.linestrings([[[0, 0], [0, 1]]] * 3),
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
    "calculation_type,breach_ids",
    [
        (CalculationType.ISOLATED, [-9999, -9999]),
        (CalculationType.CONNECTED, [1, -9999]),
        (CalculationType.DOUBLE_CONNECTED, [1, 2]),
    ],
)
def test_match_breach_ids_with_calculation_types_ok(calculation_type, breach_ids):
    nodes = Nodes(
        id=[1],
        calculation_type=calculation_type,
        breach_ids=[breach_ids],
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=[11],
    )
    nodes.breach_ids[0, :] = breach_ids

    PotentialBreachPoints.match_breach_ids_with_calculation_types(nodes)


@pytest.mark.parametrize(
    "calculation_type,breach_ids,msg",
    [
        (
            CalculationType.ISOLATED,
            [1, -9999],
            "The following objects have potential breaches, but are not (double) connected: channels [11].",
        ),
        (
            CalculationType.ISOLATED,
            [1, 2],
            "The following objects have potential breaches, but are not (double) connected: channels [11].",
        ),
        (
            CalculationType.CONNECTED,
            [1, 2],
            "The following objects have two potential breaches at the same position, but are not double connected: channels [11].",
        ),
    ],
)
def test_match_breach_ids_with_calculation_types_err(calculation_type, breach_ids, msg):
    nodes = Nodes(
        id=[1],
        calculation_type=calculation_type,
        breach_ids=[breach_ids],
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=[11],
    )
    nodes.breach_ids[0, :] = breach_ids

    with pytest.raises(SchematisationError) as e:
        PotentialBreachPoints.match_breach_ids_with_calculation_types(nodes)

    assert str(e.value) == msg
