from threedigrid_builder.base import Breaches
from threedigrid_builder.base import Levees
from threedigrid_builder.base import Lines
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import Material
from threedigrid_builder.grid import ConnectedPoints

import numpy as np
import pygeos
import pytest
from numpy.testing import assert_almost_equal, assert_equal


@pytest.fixture
def levees():
    return Levees(
        id=[1, 2],
        the_geom=[
            pygeos.linestrings([[0, 0], [10, 0], [10, 10]]),
            pygeos.linestrings([[5, 5], [8, 5]]),
        ],
        max_breach_depth=[4.0, np.nan],
        material=[Material.CLAY, -9999],
    )


@pytest.fixture
def connected_points():
    return ConnectedPoints(
        id=[0, 1, 2],
        levee_id=[1, 1, 2],
        exchange_level=[2.0, np.nan, np.nan],
    )


@pytest.fixture
def lines():
    CH = ContentType.TYPE_V2_CHANNEL
    CP = ContentType.TYPE_V2_ADDED_CALCULATION_POINT
    return Lines(
        id=[0, 1, 2, 3, 4],
        line_geometries=[
            None,
            None,
            pygeos.linestrings([[5, 5], [15, 5]]),  # crosses levee 1 at [10, 5]
            pygeos.linestrings([[3, 2], [4, 1], [5, 2]]),  # (closest on levee: [4, 0])
            pygeos.linestrings([[6, 5], [10, 5]]),  # tangent to levee 2
        ],
        content_type=[CH, -9999, CP, CP, CP],
        content_pk=[1, -9999, 0, 1, 2],
    )


def test_get_breaches(connected_points, lines, levees):
    breaches = connected_points.get_breaches(lines, levees)

    assert len(breaches) == 3

    assert_equal(breaches.id, [1, 2, 3])
    assert_equal(breaches.content_pk, [0, 1, 2])
    assert_almost_equal(breaches.coordinates, [[10, 5], [4, 0], [6, 5]])
    assert_equal(breaches.levee_id, [1, 1, 2])
    assert_equal(breaches.levl, [2, 3, 4])
    assert_almost_equal(breaches.levbr, [4, 4, np.nan])
    assert_equal(breaches.levmat, [Material.CLAY, Material.CLAY, -9999])
