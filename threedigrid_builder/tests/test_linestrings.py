import numpy as np
import pytest
import shapely
from numpy.testing import assert_almost_equal, assert_array_equal
from shapely.testing import assert_geometries_equal

from threedigrid_builder.base import Array, LinesOnLine, LineStrings, PointsOnLine


class LinearObject:
    id: int
    the_geom: shapely.Geometry


class LinearObjects(Array[LinearObject]):
    pass


@pytest.fixture
def linestrings():
    return LineStrings(
        LinearObjects(id=[0], the_geom=shapely.linestrings([[(0, 0), (6, 0), (6, 6)]]))
    )


@pytest.fixture
def two_linestrings():
    return LineStrings(
        LinearObjects(
            id=[1, 3],
            the_geom=np.array(
                [
                    shapely.linestrings([(0, 10), (10, 10)]),
                    shapely.linestrings([(0, 0), (6, 0), (6, 6)]),
                ]
            ),
        )
    )


@pytest.mark.parametrize(
    "s1d,expected_idx,expected_idx_s1d_start, expected_s1d_end",
    [
        ([], [0], [0.0], [12.0]),
        ([6.0], [0, 0], [0.0, 6.0], [6.0, 12.0]),
        ([2.0, 8.0], [0, 0, 0], [0.0, 2.0, 8.0], [2.0, 8.0, 12.0]),
    ],
)
def test_segmentize_one(
    linestrings, s1d, expected_idx, expected_idx_s1d_start, expected_s1d_end
):
    points_on_line = PointsOnLine(
        linestrings=linestrings,
        id=range(len(s1d)),
        s1d=s1d,
        linestring_idx=[0] * len(s1d),
    )
    lines = linestrings.segmentize(points_on_line)

    assert isinstance(lines, LinesOnLine)
    assert_array_equal(lines.linestring_idx, expected_idx)
    assert_array_equal(lines.s1d_start, expected_idx_s1d_start)
    assert_array_equal(lines.s1d_end, expected_s1d_end)


def test_segmentize_multiple(two_linestrings):
    points_on_line = PointsOnLine(
        linestrings=two_linestrings,
        id=[1, 2, 3],
        s1d=[5.0, 2.0, 8.0],
        linestring_idx=[0, 1, 1],
    )
    lines = two_linestrings.segmentize(points_on_line)

    assert isinstance(lines, LinesOnLine)
    assert_array_equal(lines.linestring_idx, [0, 0, 1, 1, 1])
    assert_array_equal(lines.s1d_start, [0.0, 5.0, 0.0, 2.0, 8.0])
    assert_array_equal(lines.s1d_end, [5.0, 10.0, 2.0, 8.0, 12.0])


@pytest.mark.parametrize(
    "geom,expected",
    [
        (None, None),
        ([(0, 21), (3, 42)], [(0, 21), (3, 42)]),
        ([(0, 21), (1, 2), (3, 42)], [(0, 21), (1, 2), (3, 42)]),
        ([(0, 21), (3, 42)][::-1], [(0, 21), (3, 42)]),
        ([(0, 21), (1, 2), (3, 42)][::-1], [(0, 21), (1, 2), (3, 42)]),
        ([(1, 21), (3, 41)], [(0, 21), (1, 21), (3, 41), (3, 42)]),
        ([(3, 41), (1, 21)], [(0, 21), (1, 21), (3, 41), (3, 42)]),
    ],
)
def test_sanitize(geom, expected):
    linestrings = LineStrings(
        LinearObjects(id=[0], the_geom=[shapely.linestrings(geom) if geom else None])
    )
    linestrings.sanitize(
        points_1=shapely.points([(0, 21)]), points_2=shapely.points([(3, 42)])
    )

    expected = [None] if expected is None else [shapely.linestrings(expected)]
    assert_geometries_equal(linestrings.the_geom, expected)


def test_point_on_line_the_geom_multiple(two_linestrings):
    points_on_line = PointsOnLine(
        linestrings=two_linestrings, id=[1, 2], s1d=[3.0, 4.0], linestring_idx=[0, 1]
    )
    actual = points_on_line.the_geom

    coords, index = shapely.get_coordinates(actual, return_index=True)

    assert_almost_equal(coords, [(3.0, 10.0), (4.0, 0)])
    assert_almost_equal(index, [0, 1])


def test_point_on_line_merge_with(linestrings):
    a = PointsOnLine(
        linestrings=linestrings,
        id=[0, 1],
        content_pk=[1, 2],
        s1d=[3.0, 4.0],
        linestring_idx=[0, 1],
    )
    b = PointsOnLine(
        linestrings=linestrings,
        id=[0, 1],
        content_pk=[5, 6],
        s1d=[1.0, 2.0],
        linestring_idx=[0, 1],
    )
    actual = a.merge_with(b)

    assert_array_equal(actual.s1d, [1.0, 3.0, 2.0, 4.0])
    assert_array_equal(actual.content_pk, [5, 1, 6, 2])
    assert_array_equal(actual.linestring_idx, [0, 0, 1, 1])


def test_point_on_line_at_start_at_end(two_linestrings):
    points_on_line = PointsOnLine(
        linestrings=two_linestrings,
        id=range(4),
        s1d=[0.0, 10.0, 0.0, 12.0],
        linestring_idx=[0, 0, 1, 1],
    )

    assert_array_equal(points_on_line.at_start, [True, False, True, False])
    assert_array_equal(points_on_line.at_end, [False, True, False, True])


def test_point_on_line_linestring_id(two_linestrings):
    points_on_line = PointsOnLine(
        linestrings=two_linestrings,
        id=range(3),
        linestring_idx=[0, 1, 0],
    )

    assert_array_equal(points_on_line.linestring_id, [1, 3, 1])


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
def test_line_on_line_the_geom(linestrings, start, end, expected_coords):
    lines_on_line = LinesOnLine(
        linestrings=linestrings,
        id=[1],
        s1d_start=[start],
        s1d_end=[end],
        linestring_idx=[0],
    )
    geometries = lines_on_line.the_geom
    assert geometries.shape == (1,)
    assert_almost_equal(expected_coords, shapely.get_coordinates(geometries), decimal=7)


def test_line_on_line_the_geom_multiple(two_linestrings):
    lines_on_line = LinesOnLine(
        linestrings=two_linestrings,
        id=[1, 2],
        s1d_start=[3.0, 4.0],
        s1d_end=[9.0, 10.0],
        linestring_idx=[0, 1],
    )
    segments = lines_on_line.the_geom

    coords, index = shapely.get_coordinates(segments, return_index=True)

    assert_almost_equal(
        coords, [(3.0, 10.0), (9.0, 10.0), (4.0, 0), (6.0, 0.0), (6.0, 4.0)]
    )
    assert_almost_equal(index, [0, 0, 1, 1, 1])


@pytest.mark.parametrize(
    "segment_size,expected_size,expected_count",
    [
        (3.0, 12 / 4, 3),
        (4.0, 12 / 3, 2),
        (5.0, 12 / 2, 1),
        (6.0, 12 / 2, 1),
        (8.0, 12 / 2, 1),
        (9.0, 12, 0),
        (100.0, 12, 0),
    ],
)
def test_line_on_line_interpolate_points_one(
    linestrings, segment_size, expected_size, expected_count
):
    lines_on_line = LinesOnLine(
        linestrings=linestrings,
        id=[1],
        s1d_start=[0.0],
        s1d_end=linestrings.length,
        linestring_idx=[0],
    )
    points = lines_on_line.interpolate_points(segment_size)

    assert isinstance(points, PointsOnLine)
    assert len(points) == expected_count
    assert_array_equal(points.linestring_idx, 0)
    assert_almost_equal(points.s1d, np.arange(1, expected_count + 1) * expected_size)


@pytest.mark.parametrize(
    "s1d_start,s1d_end,expected_s1d",
    [
        (0.0, 6.0, [3.0]),
        (6.0, 12.0, [9.0]),
        (1.0, 10.0, [4.0, 7.0]),
    ],
)
def test_line_on_line_interpolate_points_partial(
    linestrings, s1d_start, s1d_end, expected_s1d
):
    lines_on_line = LinesOnLine(
        linestrings=linestrings,
        id=[1],
        s1d_start=[s1d_start],
        s1d_end=[s1d_end],
        linestring_idx=[0],
    )

    points = lines_on_line.interpolate_points(3.0)

    assert isinstance(points, PointsOnLine)
    assert_array_equal(points.linestring_idx, 0)
    assert_almost_equal(points.s1d, expected_s1d)


def test_line_on_line_interpolate_points_two(two_linestrings):
    lines_on_line = LinesOnLine(
        linestrings=two_linestrings,
        id=[1, 2],
        s1d_start=[0.0, 0.0],
        s1d_end=two_linestrings.length,
        linestring_idx=[0, 1],
    )

    points = lines_on_line.interpolate_points(5.0)

    assert isinstance(points, PointsOnLine)
    assert_array_equal(points.linestring_idx, [0, 1])
    assert_almost_equal(points.s1d, [5.0, 6.0])


def test_line_on_line_interpolate_points_two_partial(linestrings):
    lines_on_line = LinesOnLine(
        linestrings=linestrings,
        id=[1, 2],
        s1d_start=[0.0, 5.0],
        s1d_end=[5.0, 12.0],
        linestring_idx=[0, 0],
    )

    points = lines_on_line.interpolate_points(3.0)

    assert isinstance(points, PointsOnLine)
    assert_array_equal(points.linestring_idx, [0, 0])
    assert_almost_equal(points.s1d, [2.5, 8.5])


@pytest.mark.parametrize(
    "s1d_1, s1d_2, expected_1, expected_2",
    [
        ([6.0], [6.0], [0], [0]),
        ([6.0], [3.0], [0], [0]),
        ([6.0], [9.0], [0], [0]),
        ([3.0, 9.0], [0.0], [0], [1]),
        ([3.0, 9.0], [3.0], [0], [1]),
        ([3.0, 9.0], [6.0], [0], [1]),
        ([3.0, 9.0], [9.0], [0], [1]),
        ([3.0, 9.0], [12.0], [0], [1]),
        ([3.0, 6.0, 9.0], [0.0], [0], [1]),
        ([3.0, 6.0, 9.0], [3.0], [0], [1]),
        ([3.0, 6.0, 9.0], [4.0], [0], [1]),
        ([3.0, 6.0, 9.0], [6.0], [0], [1]),
        ([3.0, 6.0, 9.0], [7.0], [1], [2]),
        ([3.0, 6.0, 9.0], [9.0], [1], [2]),
        ([3.0, 6.0, 9.0], [12.0], [1], [2]),
        ([0.0, 6.0, 12.0], [0.0], [0], [1]),
        ([0.0, 6.0, 12.0], [12.0], [1], [2]),
    ],
)
def test_neighbours_one(linestrings, s1d_1, s1d_2, expected_1, expected_2):
    points_1 = PointsOnLine.from_s1d(linestrings, s1d_1, [0] * len(s1d_1))
    points_2 = PointsOnLine.from_s1d(linestrings, s1d_2, [0] * len(s1d_2))

    actual_1, actual_2 = points_1.neighbours(points_2)

    assert_array_equal(actual_1, expected_1)
    assert_array_equal(actual_2, expected_2)


@pytest.mark.parametrize(
    "s1d_1, s1d_2, expected_1, expected_2",
    [
        ([0.0, 10.0, 0.0, 12.0], [5.0, 6.0], [0, 2], [1, 3]),
        ([0.0, 10.0, 0.0, 12.0], [10.0, 0.0], [0, 2], [1, 3]),
        ([0.0, 10.0, 6.0, 12.0], [10.0, 0.0], [0, 2], [1, 3]),
        ([0.0, 5.0, 0.0, 12.0], [10.0, 0.0], [0, 2], [1, 3]),
    ],
)
def test_neighbours_multiple(two_linestrings, s1d_1, s1d_2, expected_1, expected_2):
    points_1 = PointsOnLine.from_s1d(two_linestrings, s1d_1, [0, 0, 1, 1])
    points_2 = PointsOnLine.from_s1d(two_linestrings, s1d_2, [0, 1])

    actual_1, actual_2 = points_1.neighbours(points_2)

    assert_array_equal(actual_1, expected_1)
    assert_array_equal(actual_2, expected_2)
