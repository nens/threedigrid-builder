from numpy.testing import assert_almost_equal
from numpy.testing import assert_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import Channels
from threedigrid_builder.grid import CrossSectionLocations

import numpy as np
import pygeos
import pytest


@pytest.fixture
def channels():
    return Channels(
        id=np.array([51, 52]),
        the_geom=[
            pygeos.linestrings([(0, 0), (0, 8)]),
            pygeos.linestrings([(0, 0), (30, 0)]),
        ],
        # code=np.array(["one", "two"]),
        # connection_node_start_id=np.array([21, 25]),
        # connection_node_end_id=np.array([42, 33]),
        # calculation_type=np.array([101, 102]),
    )


@pytest.fixture
def channel_lines():
    return Lines(
        id=[0, 1, 2, 3],
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=[51, 52, 52, 52],
        line=[(0, 1), (2, 4), (4, 5), (5, 3)],
        line_geometries=[
            pygeos.linestrings([(0, 0), (0, 8)]),
            pygeos.linestrings([(0, 0), (10, 0)]),
            pygeos.linestrings([(10, 0), (20, 0)]),
            pygeos.linestrings([(20, 0), (30, 0)]),
        ],
        ds1d=[8, 10, 10, 10],
    )


@pytest.fixture
def cross_section_locations():
    # this fixture is closely related to channel_grid
    return CrossSectionLocations(
        id=[6, 7, 8],
        the_geom=pygeos.points([(0, 5), (30, 0), (12, 0)]),
        channel_id=[51, 52, 52],
    )


def test_set_cross_section_weights(channels, channel_lines, cross_section_locations):
    """The test set has the following situation:

    - channel with no added nodes, the cs location at the velocity point:
      O --6-- O
    - channel with 2 added nodes, with 2 cs locations:
      O --X-- O 8-X-- O --X-- 7

    The 4 velocity points cover the following situations:
    1. Only 1 CS location
    2. 2 CS locations to its right (extrapolation)
    3. In between 2 CS locations
    4. In between 2 CS locations.
    """
    cross_section_locations.apply_to_channels(channels, channel_lines)

    assert_equal(channel_lines.cross1, [6, 8, 8, 8])
    assert_equal(channel_lines.cross2, [6, 7, 7, 7])
    assert_almost_equal(channel_lines.cross_weight, [1.0, 25 / 18, 15 / 18, 5 / 18])
