from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal
from threedigrid_builder.base import LineStrings, PointsOnLine, LinesOnLine

import numpy as np
import pygeos
import pytest


@pytest.fixture
def linestrings():
    return LineStrings(id=[1], the_geom=pygeos.linestrings([[(0, 0), (6, 0), (6, 6)]]))

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
def test_interpolate_points_one(linestrings, segment_size, expected_size, expected_nodes):
    points = linestrings.interpolate_points(segment_size)

    assert isinstance(points, PointsOnLine)
    assert_array_equal(points.linestring_idx, 0)
    assert_almost_equal(points.s1d, expected_size * (np.arange(len(expected_nodes)) + 1))


def test_interpolate_points_two(two_linestrings):
    points = two_linestrings.interpolate_points(5.0)

    assert isinstance(points, PointsOnLine)
    assert_array_equal(points.linestring_idx, [0, 1])
    assert_almost_equal(points.s1d, [5.0, 6.0])


@pytest.mark.parametrize(
    "s1d,expected_idx,expected_idx_s1d_start, expected_s1d_end",
    [
        ([], [0], [0.0], [12.0]),
        ([6.0], [0, 0], [0.0, 6.0], [6.0, 12.0]),
        ([2.0, 8.0], [0, 0, 0], [0.0, 2.0, 8.0], [2.0, 8.0, 12.0]),
    ],
)
def test_segmentize_one(linestrings, s1d,expected_idx,expected_idx_s1d_start, expected_s1d_end):
    points_on_line = PointsOnLine(id=range(len(s1d)), s1d=s1d, linestring_idx=[0] * len(s1d))
    lines = linestrings.segmentize(points_on_line)

    assert isinstance(lines, LinesOnLine)
    assert_array_equal(lines.linestring_idx, expected_idx)
    assert_array_equal(lines.s1d_start, expected_idx_s1d_start)
    assert_array_equal(lines.s1d_end, expected_s1d_end)



def test_segmentize_multiple(two_linestrings):
    points_on_line = PointsOnLine(id=[1, 2, 3], s1d=[5.0, 2.0, 8.0], linestring_idx=[0, 1, 1])
    lines = two_linestrings.segmentize(points_on_line)

    assert isinstance(lines, LinesOnLine)
    assert_array_equal(lines.linestring_idx, [0, 0, 1, 1, 1])
    assert_array_equal(lines.s1d_start, [0.0, 5.0, 0.0, 2.0, 8.0])
    assert_array_equal(lines.s1d_end, [5.0, 10.0, 2.0, 8.0, 12.0])


@pytest.mark.parametrize(
    "start,end,expected_coords",
    [
        (0.0, 12.0, [(0.0, 0.0), (6.0, 0.0), (6.0, 6.0)]),
        (1.45, 12.0, [(1.45, 0.0), (6.0, 0.0), (6.0, 6.0)]),
        (6.0, 12.0, [(6.0, 0.0), (6.0, 6.0)]),
        (6.7, 12.0, [(6.0, 0.7), (6.0, 6.0)]),
        (0.0, 11.5, [(0.0, 0.0), (6.0, 0.0), (6.0, 5.5)]),
        (0.0, 6.0, [(0.0, 0.0), (6.0, 0.0)]),
        (0.0, 3.5, [(0.0, 0.0), (3.5, 0.0)]),
        (1.0, 9.0, [(1.0, 0.0), (6.0, 0.0), (6.0, 3.0)]),
        (6.0 + 1e-7, 12.0, [(6.0 + 1e-7, 0.0), (6.0, 6.0)]),
        (6.0 + 1e-8, 12.0, [(6.0, 0.0), (6.0, 6.0)]),
        (6.0, 12.0 + 1e-15, [(6.0, 0.0), (6.0, 6.0)]),
    ],
)
def test_line_on_line_as_geometries(linestrings, start, end, expected_coords):
    lines_on_line = LinesOnLine(id=[1], s1d_start=[start], s1d_end=[end], linestring_idx=[0])
    geometries = lines_on_line.as_geometries(linestrings.the_geom)
    assert geometries.shape == (1, )
    assert_almost_equal(expected_coords, pygeos.get_coordinates(geometries), decimal=7)


def test_line_on_line_as_geometries_multiple(two_linestrings):
    lines_on_line = LinesOnLine(id=[1, 2], s1d_start=[3.0, 4.0], s1d_end=[9.0, 10.0], linestring_idx=[0, 1])
    segments = lines_on_line.as_geometries(two_linestrings.the_geom)

    coords, index = pygeos.get_coordinates(segments, return_index=True)

    assert_almost_equal(
        coords, [(3.0, 10.0), (9.0, 10.0), (4.0, 0), (6.0, 0.0), (6.0, 4.0)]
    )
    assert_almost_equal(index, [0, 0, 1, 1, 1])
