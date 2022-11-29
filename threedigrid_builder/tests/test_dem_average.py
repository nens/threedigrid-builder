import numpy as np
import pytest
import shapely
from numpy.testing import assert_almost_equal

from threedigrid_builder.base import Lines, Nodes
from threedigrid_builder.grid import DemAverageAreas
from threedigrid_builder.grid.grid import Grid


@pytest.fixture
def grid():
    id = np.array([0, 1, 2, 3])
    calculation_type = np.array([-9999, -9999, -9999, -9999])
    bounds = np.array(
        [
            [0.0, 0.0, 1.0, 1.0],
            [0.0, 1.0, 1.0, 2.0],
            [1.0, 0.0, 2.0, 1.0],
            [1.0, 1.0, 2.0, 2.0],
        ],
        dtype=np.float64,
    )
    nodes = Nodes(id=id, calculation_type=calculation_type, bounds=bounds)
    grid = Grid(nodes=nodes, lines=Lines(id=np.array([1])))
    grid._cell_tree = shapely.STRtree(shapely.box(*grid.nodes.bounds.T))
    return grid


@pytest.fixture
def dem_average_areas():
    return DemAverageAreas(
        id=np.array([1], dtype=np.int32),
        the_geom=np.array([shapely.box(0.2, 0.2, 0.8, 1.8)]),
    )


@pytest.fixture
def dem_average_areas_no_intersect():
    return DemAverageAreas(
        id=np.array([1], dtype=np.int32),
        the_geom=np.array([shapely.box(11, 11, 18, 18)]),
    )


def test_dem_average_areas(grid, dem_average_areas):
    grid.set_dem_averaged_cells(dem_average_areas)
    assert_almost_equal(grid.nodes.has_dem_averaged, [1, 1, -9999, -9999])


def test_dem_average_areas_no_intersections(grid, dem_average_areas_no_intersect):
    grid.set_dem_averaged_cells(dem_average_areas_no_intersect)
    assert_almost_equal(grid.nodes.has_dem_averaged, [-9999, -9999, -9999, -9999])
