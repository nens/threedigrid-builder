import numpy as np
import pytest
from numpy.testing import assert_equal

from threedigrid_builder.base import Lines
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import Windshieldings


@pytest.fixture
def channel_lines():
    # ordered by [content_pk, position on channel]
    return Lines(
        id=[0, 1, 2, 3, 4, 5, 6],
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=[51, 52, 53, 53, 53, 54, 54],
        line=[(0, 1), (0, 2), (0, 3), (3, 4), (4, 5), (0, 6), (6, 7)],
    )


@pytest.fixture
def mixed_lines():
    # ordered by [content_pk, position on channel]
    return Lines(
        id=[0, 1, 2, 3, 4, 5, 6],
        content_type=[
            ContentType.TYPE_V2_PIPE,
            ContentType.TYPE_V2_CULVERT,
            ContentType.TYPE_V2_CHANNEL,
            ContentType.TYPE_V2_CHANNEL,
            ContentType.TYPE_V2_CHANNEL,
            ContentType.TYPE_V2_WEIR,
            ContentType.TYPE_V2_WEIR,
        ],
        content_pk=[53, 53, 53, 53, 53, 53, 53],
        line=[(0, 1), (0, 2), (0, 3), (3, 4), (4, 5), (0, 6), (6, 7)],
    )


@pytest.fixture
def different_lines():
    # ordered by [content_pk, position on channel]
    return Lines(
        id=[0, 1, 2, 3, 4, 5, 6],
        content_type=[
            ContentType.TYPE_V2_PIPE,
            ContentType.TYPE_V2_CULVERT,
            ContentType.TYPE_V2_CHANNEL,
            ContentType.TYPE_V2_CHANNEL,
            ContentType.TYPE_V2_CHANNEL,
            ContentType.TYPE_V2_WEIR,
            ContentType.TYPE_V2_WEIR,
        ],
        content_pk=[53, 53, 1, 2, 3, 53, 53],
        line=[(0, 1), (0, 2), (0, 3), (3, 4), (4, 5), (0, 6), (6, 7)],
    )


@pytest.fixture
def windshieldings():
    north = np.array([0.1, 0.2, 0.5, 0.4], dtype=np.float64)
    east = np.array([0.0, 0.1, 0.3, 0.9], dtype=np.float64)
    south = np.array([0.3, 0.4, 0.1, 0.8], dtype=np.float64)
    west = np.array([0.4, 0.3, 0.2, 0.1], dtype=np.float64)
    return Windshieldings(
        id=[1, 2, 3, 4],
        channel_id=[51, 53, 52, 54],
        north=north,
        northeast=north + 1,
        east=east,
        southeast=east + 1,
        south=south,
        southwest=south + 1,
        west=west,
        northwest=west + 1,
    )


@pytest.fixture
def windshieldings_channel_selection():
    north = np.array([0.1], dtype=np.float64)
    east = np.array([0.0], dtype=np.float64)
    south = np.array([0.3], dtype=np.float64)
    west = np.array([0.4], dtype=np.float64)
    return Windshieldings(
        id=[1],
        channel_id=[51],
        north=north,
        northeast=north + 1,
        east=east,
        southeast=east + 1,
        south=south,
        southwest=south + 1,
        west=west,
        northwest=west + 1,
    )


def test_set_windshielding(channel_lines, windshieldings):
    windshieldings.apply_to_lines(channel_lines)
    expected_windshieldings = np.array(
        [
            [0.1, 1.1, 0.0, 1.0, 0.3, 1.3, 0.4, 1.4],
            [0.5, 1.5, 0.3, 1.3, 0.1, 1.1, 0.2, 1.2],
            [0.2, 1.2, 0.1, 1.1, 0.4, 1.4, 0.3, 1.3],
            [0.2, 1.2, 0.1, 1.1, 0.4, 1.4, 0.3, 1.3],
            [0.2, 1.2, 0.1, 1.1, 0.4, 1.4, 0.3, 1.3],
            [0.4, 1.4, 0.9, 1.9, 0.8, 1.8, 0.1, 1.1],
            [0.4, 1.4, 0.9, 1.9, 0.8, 1.8, 0.1, 1.1],
        ],
        dtype=np.float64,
    )
    assert_equal(channel_lines.windshieldings, expected_windshieldings)


def test_set_windshielding_no_windshieldings(channel_lines, windshieldings):
    Windshieldings(id=[]).apply_to_lines(channel_lines)
    expected_windshieldings = np.full((7, 8), np.nan, dtype=np.float64)
    assert_equal(channel_lines.windshieldings, expected_windshieldings)


def test_set_windshielding_different_lines(mixed_lines, windshieldings):
    windshieldings.apply_to_lines(mixed_lines)
    expected_windshieldings = np.array(
        [
            [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            [0.2, 1.2, 0.1, 1.1, 0.4, 1.4, 0.3, 1.3],
            [0.2, 1.2, 0.1, 1.1, 0.4, 1.4, 0.3, 1.3],
            [0.2, 1.2, 0.1, 1.1, 0.4, 1.4, 0.3, 1.3],
            [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
        ],
        dtype=np.float64,
    )
    assert_equal(mixed_lines.windshieldings, expected_windshieldings)


def test_set_windshielding_channel_selection(
    channel_lines, windshieldings_channel_selection
):
    windshieldings_channel_selection.apply_to_lines(channel_lines)
    expected_windshieldings = np.array(
        [
            [0.1, 1.1, 0.0, 1.0, 0.3, 1.3, 0.4, 1.4],
            [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
            [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan],
        ],
        dtype=np.float64,
    )
    assert_equal(channel_lines.windshieldings, expected_windshieldings)
