from numpy.testing import assert_almost_equal
from numpy.testing import assert_equal
from threedigrid_builder.base import Levees
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import LineType
from threedigrid_builder.grid import Obstacles
from threedigrid_builder.grid.grid import Grid

import numpy as np
import pygeos
import pytest


@pytest.fixture
def grid():
    id = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], dtype=np.int32)
    kcu = np.array([98, 98, 98, 98, 98, 99, 99, 99, 99, 99], dtype=np.int32)
    line_coords = np.array(
        [
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
        dtype=np.float64,
    )
    lines = Lines(id=id, kcu=kcu, line_coords=line_coords)
    return Grid(nodes=Nodes(id=np.array([1])), lines=lines)


@pytest.fixture
def obstacles():
    return Obstacles(
        id=np.array([1, 2], dtype=np.int32),
        crest_level=np.array([-0.1, -0.3], dtype=np.float64),
        the_geom=pygeos.linestrings(
            [[[12.0, 11.0], [17.5, 17.0]], [[13.0, 11.0], [13.0, 17.0]]]
        ),
    )


@pytest.fixture
def obstacles_no_intersect():
    return Obstacles(
        id=[1],
        crest_level=[-0.1],
        the_geom=pygeos.linestrings(
            [[[11, 11], [11, 18], [18, 18], [18, 11], [11, 11]]]
        ),
    )


def test_set_obstacles(grid, obstacles):
    grid.set_obstacles(obstacles, Levees(id=[]))

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
    grid.set_obstacles(obstacles_no_intersect, Levees(id=[]))

    assert np.all(grid.lines.kcu != LineType.LINE_2D_OBSTACLE)
    assert_equal(grid.lines.flod, np.nan)
    assert_equal(grid.lines.flou, np.nan)
