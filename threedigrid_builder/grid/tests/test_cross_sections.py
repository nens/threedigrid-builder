from numpy.testing import assert_almost_equal
from numpy.testing import assert_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import Channels
from threedigrid_builder.grid import CrossSectionLocations
from threedigrid_builder.grid.cross_sections import compute_weights

import pygeos
import pytest


@pytest.fixture
def channels():
    return Channels(
        id=[51, 52, 53, 54],
        the_geom=[
            pygeos.linestrings([(0, 0), (3, 0)]),  # 1 segment of size 3
            pygeos.linestrings([(55, 3), (60, 3)]),  # 1 segment of size 5
            pygeos.linestrings([(3, 1), (36, 1)]),  # 3 segments of size 11
            pygeos.linestrings([(40, 2), (54, 2)]),  # 2 segments of size 7
        ],
    )


@pytest.fixture
def channel_lines():
    # ordered by [content_pk, position on channel]
    return Lines(
        id=[0, 1, 2, 3, 4, 5, 6],
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=[51, 52, 53, 53, 53, 54, 54],
        ds1d=[3, 5, 11, 11, 11, 7, 7],
    )


@pytest.fixture
def cross_section_locations():
    return CrossSectionLocations(
        id=[1, 2, 3, 5, 6],
        the_geom=pygeos.points([(1, 0), (26, 1), (16, 1), (59, 3), (49, 2)]),
        channel_id=[51, 53, 53, 52, 54],
    )


def test_compute_weights(cross_section_locations, channels, channel_lines):
    """The test set has the following situation:

    - channel with no added nodes and with 1 cs location
      O 1-X- O
    - channel with no added nodes and with 2 cs locations:
      O --X-5O
    - channel with 2 added nodes, with 2 cs locations, extrapolation:
      O --X-- O 3-X-- O 2-X-- O
    - channel with 1 added node, with 1 cs location, no extrapolation:
      O --X-- O 6-X-- O
    Here O means connection nodes, --X-- a line with a velocity point, and
    CS locations are given as numbers (ids).

    The 7 lines (velocity points) cover all possible situations:
    0. Only 1 CS location in the channel, extrapolation would lead to index -1
    1. Extrapolation on the left
    2. Interpolation
    3. Extrapolation on the right
    4. Equalization on the left
    5. Equalization on the right
    6. Only 1 CS location in the channel, extrapol. lead to too large index
    """
    expected_cross1 = [1, 5, 3, 3, 3, 6, 6]
    expected_cross2 = [1, 5, 2, 2, 2, 6, 6]
    expected_weight = [1.0, 1.0, 1.75, 0.65, -0.45, 1.0, 1.0]

    cross1, cross2, cross_weights = compute_weights(
        channel_lines.content_pk, channel_lines.ds1d, cross_section_locations, channels
    )

    assert_equal(cross1, expected_cross1)
    assert_equal(cross2, expected_cross2)
    assert_almost_equal(cross_weights, expected_weight)
