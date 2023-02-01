import numpy as np
import pytest
import shapely
from numpy.testing import assert_array_equal

from threedigrid_builder.grid import Channels, CrossSectionLocations


@pytest.fixture
def channels():
    return Channels(
        the_geom=shapely.linestrings([[(0, 0), (6, 0), (6, 6)]]),
        dist_calc_points=np.array([5.0]),
        id=np.array([1]),
        code=np.array(["one"]),
        connection_node_start_id=np.array([21]),
        connection_node_end_id=np.array([42]),
        calculation_type=np.array([2]),
    )


def test_get_1d2d_exchange_levels(channels):
    locations = CrossSectionLocations(
        id=[2, 5],
        the_geom=shapely.points([(0, 0), [6, 6]]),
        bank_level=[1.0, 13.0],
        channel_id=[1, 1],
    )

    actual = channels.get_1d2d_exchange_levels(
        content_pk=np.array([1, 1]), s1d=np.array([3.0, 6.0]), locations=locations
    )

    # bank levels are interpolated
    assert_array_equal(actual, [4.0, 7.0])


def test_is_closed(channels):
    actual = channels.is_closed(np.array([1, 3]))

    assert_array_equal(actual, [False, False])


@pytest.mark.parametrize(
    "thickness,hc_out,hc_in,expected",
    [
        (0.1, 3.0, 2.0, True),
        (0.0, 3.0, 2.0, False),
        (np.nan, 3.0, 2.0, False),
        (0.1, np.nan, 2.0, False),
        (0.1, 3.0, np.nan, False),
    ],
)
def test_has_groundwater_exchange(thickness, hc_out, hc_in, expected):
    channels = Channels(
        id=[1],
        exchange_thickness=thickness,
        hydraulic_conductivity_out=hc_out,
        hydraulic_conductivity_in=hc_in,
    )

    actual = channels.has_groundwater_exchange

    assert len(actual) == 1
    assert actual[0] == expected
