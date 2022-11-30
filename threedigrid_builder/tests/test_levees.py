import numpy as np
import pygeos
import pytest
from numpy.testing import assert_almost_equal, assert_equal

from threedigrid_builder.base import Lines, LineStrings, PointsOnLine
from threedigrid_builder.constants import ContentType, Material
from threedigrid_builder.grid import (
    Breaches,
    ConnectedPoints,
    Levees,
    PotentialBreaches,
)


@pytest.fixture
def levees():
    return Levees(
        id=[1, 2, 3],
        the_geom=[
            pygeos.linestrings([[0, 0], [10, 0], [10, 10]]),
            pygeos.linestrings([[5, 5], [8, 5]]),
            pygeos.linestrings([[0, 0], [1, 1]]),
        ],
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
    CH = ContentType.TYPE_V2_CHANNEL
    CP = ContentType.TYPE_V2_ADDED_CALCULATION_POINT
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

    assert isinstance(breaches, Breaches)
    assert len(breaches) == 3

    assert_equal(breaches.id, [0, 1, 2])
    assert_equal(breaches.content_pk, [0, 1, 2])
    assert_almost_equal(breaches.coordinates, [[10, 5], [4, 0], [6, 5]])
    assert_equal(breaches.levee_id, [1, 1, 2])
    assert_equal(breaches.levl, [2, 3, 4])
    assert_almost_equal(breaches.levbr, [4, 4, 2])
    assert_equal(breaches.levmat, [Material.CLAY, Material.CLAY, Material.SAND])


def test_no_breaches(connected_points, lines, levees):
    connected_points.levee_id[:] = -9999
    breaches = connected_points.get_breaches(lines, levees)

    assert isinstance(breaches, Breaches)
    assert len(breaches) == 0


def test_potential_breach_sides():
    potential_breaches = PotentialBreaches(
        id=[1],
        channel_id=[1],
        the_geom=[pygeos.linestrings([[0, 0], [0, 1]])],
    )
    assert_almost_equal(pygeos.get_coordinates(potential_breaches.side_1d), [[0, 0]])
    assert_almost_equal(pygeos.get_coordinates(potential_breaches.side_2d), [[0, 1]])


def test_potential_breach_merge():
    before = PointsOnLine(
        linestrings=LineStrings(pygeos.linestrings([[[0, 0], [10, 0]]])),
        id=[0, 1, 2, 3, 4, 5],
        s1d=[0.0, 2.0, 5.0, 6.0, 7.0, 10.0],
        linestring_idx=[0, 0, 0, 0, 0, 0],
    )
    after = PotentialBreaches.merge(before, tolerance=2.0)

    assert_almost_equal(after.s1d, [0.0, 2.0, 5.5, 10.0])
    assert_almost_equal(after.linestring_idx, [0, 0, 0, 0])
