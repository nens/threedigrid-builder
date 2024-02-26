import numpy as np
import pytest
import shapely
from numpy.testing import assert_almost_equal, assert_equal

from threedigrid_builder.base import Lines, Nodes
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import (
    Channels,
    CrossSectionDefinitions,
    CrossSectionLocations,
)
from threedigrid_builder.grid.cross_section_locations import (
    compute_bottom_level,
    compute_weights,
)


@pytest.fixture
def channels():
    return Channels(
        id=[51, 52, 53, 54],
        the_geom=[
            shapely.linestrings([(0, 0), (3, 0)]),  # 1 segment of size 3
            shapely.linestrings([(55, 3), (60, 3)]),  # 1 segment of size 5
            shapely.linestrings([(3, 1), (36, 1)]),  # 3 segments of size 11
            shapely.linestrings([(40, 2), (54, 2)]),  # 2 segments of size 7
        ],
    )


@pytest.fixture
def channel_nodes():
    CH = ContentType.TYPE_V2_CHANNEL
    CN = ContentType.TYPE_V2_CONNECTION_NODES
    return Nodes(
        id=range(8),
        content_type=[CN, CN, CN, CH, CH, CN, CH, CN],
        dmax=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
    )


@pytest.fixture
def channel_lines():
    # ordered by [content_pk, position on channel]
    return Lines(
        id=[0, 1, 2, 3, 4, 5, 6],
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=[51, 52, 53, 53, 53, 54, 54],
        ds1d=[3, 5, 11, 11, 11, 7, 7],
        s1d=[1.5, 2.5, 5.5, 16.5, 27.5, 3.5, 10.5],
        line=[(0, 1), (0, 2), (0, 3), (3, 4), (4, 5), (0, 6), (6, 7)],
        dpumax=[3.0, 3.0, 0, 0, 0, 0, 0],
    )


@pytest.fixture
def locations():
    return CrossSectionLocations(
        id=[1, 2, 3, 5, 6],
        the_geom=shapely.points([(1, 0), (26, 1), (16, 1), (58, 3), (49, 2)]),
        channel_id=[51, 53, 53, 52, 54],
        friction_type=[1, 1, 1, 2, 2],
        friction_value=[30, 35, 40, 0.02, 0.03],
        vegetation_stem_density=[1, 0.3, 0.4, np.nan, np.nan],
        vegetation_stem_diameter=[0.5, 0.5, 0.5, np.nan, np.nan],
        vegetation_height=[0.2, 1, 1, np.nan, np.nan],
        vegetation_drag_coefficient=[1, 1, 1, np.nan, np.nan],
        reference_level=[1.0, 2.0, 3.0, 5.0, 6.0],
        definition_id=[3, 3, 4, 4, 5],
    )


@pytest.fixture
def definitions():
    return CrossSectionDefinitions(id=[3, 4, 5])


@pytest.mark.parametrize("extrapolate", [True, False])
def test_compute_weights(locations, channels, channel_lines, extrapolate):
    """The test set has the following situation:

    - channel 51 with no added nodes and with 1 cs location
      O 1-X-- O
    - channel 52 with no added nodes and with 1 cs location:
      O --X-5 O
    - channel 53 with 2 added nodes, with 2 cs locations, extrapolation:
      O --X-- O 2-X-- O 3-X-- O
    - channel 54 with 1 added node, with 1 cs location, no extrapolation:
      O --X-- O 6-X-- O
    Here O means connection nodes, --X-- a line with a velocity point, and
    CS locations are given as numbers (ids).

    The 7 lines (velocity points) cover all possible situations:
    0. Only 1 CS location in the channel, extrapolation would lead to index -1
    1. Equalization on the right
    2. Extrapolation on the left
    3. Interpolation
    4. Extrapolation on the right
    5. Equalization on the left
    6. Only 1 CS location in the channel, extrapol. lead to too large index
    """
    expected_cross_loc1 = [1, 5, 3, 3, 3, 6, 6]
    expected_cross_loc2 = [1, 5, 2, 2, 2, 6, 6]
    expected_weight = [1.0, 1.0, 1.75, 0.65, -0.45, 1.0, 1.0]
    if not extrapolate:
        expected_weight = np.clip(expected_weight, 0, 1)

    cross_loc1, cross_loc2, cross_weights = compute_weights(
        channel_lines.content_pk, channel_lines.s1d, locations, channels, extrapolate
    )

    assert_equal(cross_loc1, expected_cross_loc1)
    assert_equal(cross_loc2, expected_cross_loc2)
    assert_almost_equal(cross_weights, expected_weight)


def test_compute_weights_edge_effects(channels):
    locations = CrossSectionLocations(
        id=range(3),
        the_geom=shapely.points([(0, 0), (55, 3), (3, 0)]),
        channel_id=[51, 52, 51],
    )
    expected_cross_loc1 = [0, 0, 1, 1]
    expected_cross_loc2 = [2, 2, 1, 1]
    expected_weight = [1.0, 0.0, 1.0, 1.0]

    cross_loc1, cross_loc2, cross_weights = compute_weights(
        [51, 51, 52, 52],
        [0.0, 3.0, 0.0, 55.0],
        locations,
        channels[:2],
        extrapolate=True,
    )

    assert_equal(cross_loc1, expected_cross_loc1)
    assert_equal(cross_loc2, expected_cross_loc2)
    assert_almost_equal(cross_weights, expected_weight)


@pytest.mark.parametrize(
    "extrapolate,expected",
    [
        (True, [1.0, 5.0, 3.75, 2.65, 1.55, 6.0, 6.0]),
        (False, [1.0, 5.0, 3.0, 2.65, 2.0, 6.0, 6.0]),
    ],
)
def test_compute_bottom_level(
    locations, channels, channel_lines, extrapolate, expected
):
    """Same setup as test_compute_weights, but now testing the derived bottom levels"""
    actual = compute_bottom_level(
        channel_lines.content_pk,
        channel_lines.s1d,
        locations,
        channels,
        extrapolate=extrapolate,
    )

    assert_almost_equal(actual, expected)


def test_apply_to_lines(channels, channel_lines, locations):
    locations.apply_to_lines(channel_lines, channels, extrapolate=False)

    assert_equal(channel_lines.cross_id1, [3, 4, 4, 4, 4, 5, 5])
    assert_equal(channel_lines.cross_id2, [3, 4, 3, 3, 3, 5, 5])
    assert_equal(channel_lines.frict_type1, [1, 2, 1, 1, 1, 2, 2])
    assert_equal(channel_lines.frict_type2, [1, 2, 1, 1, 1, 2, 2])
    assert_equal(channel_lines.frict_value1, [30.0, 0.02, 40.0, 40.0, 40.0, 0.03, 0.03])
    assert_equal(channel_lines.frict_value2, [30.0, 0.02, 35.0, 35.0, 35.0, 0.03, 0.03])
    assert_almost_equal(channel_lines.veg_coef1, [0.5, 0.0, 0.2, 0.2, 0.2, 0.0, 0.0])
    assert_almost_equal(channel_lines.veg_coef2, [0.5, 0.0, 0.15, 0.15, 0.15, 0.0, 0.0])
    assert_almost_equal(channel_lines.veg_height1, [0.2, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0])
    assert_almost_equal(channel_lines.veg_height2, [0.2, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0])
    assert_almost_equal(
        channel_lines.cross_weight, [1.0, 1.0, 1.0, 0.65, 0.0, 1.0, 1.0]
    )
    assert_almost_equal(
        channel_lines.invert_level_start_point, [1.0, 5.0, 3.0, 3.0, 2.1, 6.0, 6.0]
    )
    assert_almost_equal(
        channel_lines.invert_level_end_point, [1.0, 5.0, 3.0, 2.1, 2.0, 6.0, 6.0]
    )
    assert_almost_equal(channel_lines.dpumax, [1.0, 5.0, 3, 3.0, 2.1, 6.0, 6.0])


def test_apply_to_lines_extrapolate(channels, channel_lines, locations):
    locations.apply_to_lines(channel_lines, channels, extrapolate=True)

    assert_almost_equal(
        channel_lines.cross_weight, [1.0, 1.0, 1.75, 0.65, -0.45, 1.0, 1.0]
    )
    assert_almost_equal(
        channel_lines.invert_level_start_point, [1.0, 5.0, 4.3, 3.2, 2.1, 6.0, 6.0]
    )
    assert_almost_equal(
        channel_lines.invert_level_end_point, [1.0, 5.0, 3.2, 2.1, 1.0, 6.0, 6.0]
    )
    assert_almost_equal(channel_lines.dpumax, [1.0, 5.0, 4.3, 3.2, 2.1, 6.0, 6.0])


def test_apply_to_lines_project_locations(channels, channel_lines, locations):
    # CrossSectionLocations are not on the channels
    # - channel 51 (cs index 0): displaced 1 (should be projected on the line)
    # - channel 52 (cs index 3): displaced far away (should be projected to line end)
    locations.the_geom[0] = shapely.points((1, 1))  # id = 1
    locations.the_geom[3] = shapely.points((1000, 3))  # id = 5
    locations.apply_to_lines(channel_lines, channels, extrapolate=False)

    assert_equal(channel_lines.cross_id1, [3, 4, 4, 4, 4, 5, 5])
    assert_equal(channel_lines.cross_id2, [3, 4, 3, 3, 3, 5, 5])
    assert_almost_equal(
        channel_lines.cross_weight, [1.0, 1.0, 1.0, 0.65, 0.0, 1.0, 1.0]
    )
