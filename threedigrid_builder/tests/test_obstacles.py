import numpy as np
import pytest
import shapely
from numpy.testing import assert_almost_equal, assert_equal

from threedigrid_builder.base import Lines, Nodes
from threedigrid_builder.constants import LineType
from threedigrid_builder.grid import Obstacles
from threedigrid_builder.grid.grid import Grid


@pytest.fixture
def lines():
    lines = Lines(
        id=range(10),
        kcu=[98, 98, 98, 98, 98, 99, 99, 99, 99, 99],
        line_coords=[
            [12.0, 12.0, 16.0, 12.0],
            [12.0, 16.0, 15.0, 15.0],
            [12.0, 16.0, 15.0, 17.0],
            [15.0, 15.0, 17.0, 15.0],
            [15.0, 17.0, 17.0, 17.0],
            [12.0, 12.0, 12.0, 16.0],
            [16.0, 12.0, 15.0, 15.0],
            [16.0, 12.0, 17.0, 15.0],
            [15.0, 15.0, 15.0, 17.0],
            [17.0, 15.0, 17.0, 17.0],
        ],
    )
    lines.fix_line_geometries()
    return lines


@pytest.fixture
def grid(lines):
    return Grid(nodes=Nodes(id=np.array([1])), lines=lines)


@pytest.fixture
def obstacles():
    return Obstacles(
        id=[1, 2],
        crest_level=[-0.1, -0.3],
        the_geom=shapely.linestrings(
            [[[12.0, 11.0], [17.5, 17.0]], [[13.0, 11.0], [13.0, 17.0]]]
        ),
    )


@pytest.fixture
def obstacles_no_intersect():
    return Obstacles(
        id=[1],
        crest_level=[-0.1],
        the_geom=shapely.linestrings(
            [[[11, 11], [11, 18], [18, 18], [18, 11], [11, 11]]]
        ),
    )


def test_set_obstacles(grid, obstacles):
    grid.set_obstacles(obstacles)

    expected_kcu = np.array(
        [102, 102, 102, 102, 98, 99, 103, 99, 99, 103], dtype=np.int32
    )
    assert_equal(grid.lines.kcu, expected_kcu)
    expected_flodu = np.array(
        [-0.1, -0.3, -0.3, -0.1, np.nan, np.nan, -0.1, np.nan, np.nan, -0.1],
        dtype=np.float64,
    )
    assert_almost_equal(grid.lines.flod, expected_flodu)
    assert_almost_equal(grid.lines.flou, expected_flodu)


def test_set_obstacles_no_intersect(grid, obstacles_no_intersect):
    grid.set_obstacles(obstacles_no_intersect)

    assert np.all(grid.lines.kcu != LineType.LINE_2D_OBSTACLE)
    assert_equal(grid.lines.flod, np.nan)
    assert_equal(grid.lines.flou, np.nan)


def test_compute_dpumax(lines: Lines, obstacles: Obstacles):
    expected = [-0.1, -0.3, -0.3, -0.1, np.nan, np.nan, -0.1, np.nan, np.nan, -0.1]
    actual, idx = obstacles.compute_dpumax(lines, where=np.arange(len(lines)))

    assert_equal(idx, [0, 1, 1, 0, -9999, -9999, 0, -9999, -9999, 0])
    assert_almost_equal(actual, expected)


def test_compute_dpumax_no_intersections(
    lines: Lines, obstacles_no_intersect: Obstacles
):
    actual, idx = obstacles_no_intersect.compute_dpumax(
        lines, where=np.arange(len(lines))
    )

    assert len(actual) == len(lines)
    assert_equal(idx, -9999)
    assert_almost_equal(actual, np.nan)


def test_compute_dpumax_no_obstacles(lines: Lines):
    actual, idx = Obstacles(id=[]).compute_dpumax(lines, where=np.arange(len(lines)))

    assert len(actual) == len(lines)
    assert_equal(idx, -9999)
    assert_almost_equal(actual, np.nan)


def test_compute_dpumax_where(lines: Lines, obstacles: Obstacles):
    expected = [-0.1, -0.1, np.nan, -0.1]
    actual, idx = obstacles.compute_dpumax(lines, where=[0, 3, 4, 6])

    assert_equal(idx, [0, 0, -9999, 0])
    assert_almost_equal(actual, expected)
