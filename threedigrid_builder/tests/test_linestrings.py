from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal
from threedigrid_builder.base import LineStrings, PointsOnLine

import numpy as np
import pygeos
import pytest


@pytest.fixture
def two_linestrings():
    return LineStrings(id=[1, 2], the_geom=[
            pygeos.linestrings([(0, 10), (10, 10)]),
            pygeos.linestrings([(0, 0), (6, 0), (6, 6)]),
        ]
    )


@pytest.mark.parametrize(
    "segment_size,expected_size,expected_nodes",
    [
        (3.0, 12 / 4, [(3, 0), (6, 0), (6, 3)]),
        (4.0, 12 / 3, [(4, 0), (6, 2)]),
        (5.0, 12 / 2, [(6, 0)]),
        (6.0, 12 / 2, [(6, 0)]),
        (8.0, 12 / 2, [(6, 0)]),
        (9.0, 12, np.empty((0, 2), dtype=float)),
        (100.0, 12, np.empty((0, 2), dtype=float)),
    ],
)
def test_segmentize_one(segment_size, expected_size, expected_nodes):
    line = LineStrings(id=[1], the_geom=pygeos.linestrings([[(0, 0), (6, 0), (6, 6)]]))
    points = line.segmentize(segment_size)

    assert isinstance(points, PointsOnLine)
    assert_array_equal(points.linestring_idx, 0)
    assert_almost_equal(points.s1d, expected_size * (np.arange(len(expected_nodes)) + 1))


def test_segmentize_two(two_linestrings):
    points = two_linestrings.segmentize(5.0)

    assert isinstance(points, PointsOnLine)
    assert_array_equal(points.linestring_idx, [0, 1])
    assert_almost_equal(points.s1d, [5.0, 6.0])
