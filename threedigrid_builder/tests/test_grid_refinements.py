import numpy as np
import pytest
import shapely
from numpy.testing import assert_array_equal

from threedigrid_builder.grid import GridRefinements

EPSILON = 1e-4


@pytest.fixture
def refinements():
    return GridRefinements(
        id=[5],
        refinement_level=[1],
        the_geom=[None],
    )


def test_rasterize_no_refinement():
    refinements = GridRefinements(id=[])
    actual = refinements.rasterize(
        origin=(0, 0), height=4, width=3, cell_size=0.5, no_data_value=99
    )
    assert actual.shape == (4, 3)
    assert actual.dtype == np.int32
    assert_array_equal(actual, 99)


def test_rasterize_none_refinement(refinements):
    actual = refinements.rasterize(
        origin=(0, 0), height=4, width=3, cell_size=0.5, no_data_value=99
    )
    assert actual.shape == (4, 3)
    assert actual.dtype == np.int32
    assert_array_equal(actual, 99)


def test_rasterize_outside(refinements):
    refinements.the_geom[0] = shapely.linestrings([[130.1, 11.9], [130.9, 11.9]])
    actual = refinements.rasterize(
        origin=(0, 0), height=4, width=3, cell_size=0.5, no_data_value=99
    )
    assert actual.shape == (4, 3)
    assert actual.dtype == np.int32
    assert_array_equal(actual, 99)


def test_rasterize_line_refinement(refinements):
    refinements.the_geom[0] = shapely.linestrings(
        [[15.0, 17.5], [17.5, 17.5], [17.5, 14.0]]
    )
    actual = refinements.rasterize(
        origin=(12.0, 10.0), height=8, width=6, cell_size=1.0, no_data_value=3
    )
    # Notes:
    # - origin == topleft (x, y order)
    # - geom starts at x=15.0, and does not touch cell (7, 2) which covers x=[14.0, 15.0)
    # - geom ends at y=14.0, and does not touch cell (3, 5) which covers x=[13.0, 14.0)
    expected = np.array(
        [
            [3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 1],
            [3, 3, 3, 3, 3, 1],
            [3, 3, 3, 3, 3, 1],
            [3, 3, 3, 1, 1, 1],
        ]
    )
    assert_array_equal(actual, expected)


def test_rasterize_line_refinement2(refinements):
    refinements.the_geom[0] = shapely.linestrings(
        [[15.0, 14.5], [12.5, 14.5], [12.5, 16.0]]
    )
    actual = refinements.rasterize(
        origin=(12.0, 10.0), height=8, width=6, cell_size=1.0, no_data_value=3
    )
    # Notes:
    # - origin == topleft (x, y order)
    # - geom starts at x=15.0, and does touch cell (4, 3) which covers x=[15.0, 16.0)
    # - geom ends at y=16.0, and does touch cell (6, 0) which covers y=[16.0, 17.0)
    expected = np.array(
        [
            [3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3],
            [1, 1, 1, 1, 3, 3],
            [1, 3, 3, 3, 3, 3],
            [1, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3],
        ]
    )
    assert_array_equal(actual, expected)


def test_rasterize_poly_refinement(refinements):
    refinements.the_geom[0] = shapely.box(
        13.0 + EPSILON, 11.0 + EPSILON, 14.0 + EPSILON, 13.0 + EPSILON
    )
    actual = refinements.rasterize(
        origin=(12.0, 10.0), height=8, width=6, cell_size=0.5, no_data_value=3
    )
    # Notes:
    # - origin == topleft (x, y order)
    # - cell_size == 0.5 in this test
    # - geom covers x=[13.0, 14.0], so it excludes column 1 [12.5, 13.0)
    #   and includes column 4 [14.0, 14.5)
    # - geom covers y=[11.0, 13.0] so it excludes row 1 [10.5, 11.0) and
    #   includes row 6 [13.0, 13.5)
    expected = np.array(
        [
            [3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3],
            [3, 3, 1, 1, 1, 3],
            [3, 3, 1, 1, 1, 3],
            [3, 3, 1, 1, 1, 3],
            [3, 3, 1, 1, 1, 3],
            [3, 3, 1, 1, 1, 3],
            [3, 3, 3, 3, 3, 3],
        ]
    )
    assert_array_equal(actual, expected)


def test_quadtree_small_line_refinement(refinements):
    refinements.the_geom[0] = shapely.linestrings([[13.1, 11.9], [13.9, 11.9]])
    actual = refinements.rasterize(
        origin=(12.0, 10.0), height=8, width=6, cell_size=1.0, no_data_value=3
    )
    expected = np.array(
        [
            [3, 3, 3, 3, 3, 3],
            [3, 1, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3],
            [3, 3, 3, 3, 3, 3],
        ]
    )
    assert_array_equal(actual, expected)


@pytest.mark.parametrize("reverse", [False, True])
def test_rasterize_multiple_refinements(reverse):
    refinements = GridRefinements(
        id=[5, 8],
        refinement_level=[1, 2],
        the_geom=[
            shapely.box(13.0 + EPSILON, 11.0 + EPSILON, 14.0 + EPSILON, 13.0 + EPSILON),
            shapely.box(12.5 + EPSILON, 10.5 + EPSILON, 13.5 + EPSILON, 12.5 + EPSILON),
        ],
    )
    if reverse:
        refinements.refinement_level = refinements.refinement_level[::-1]
        refinements.the_geom = refinements.the_geom[::-1]
    actual = refinements.rasterize(
        origin=(12.0, 10.0), height=8, width=6, cell_size=0.5, no_data_value=3
    )
    # Notes:
    # - see test_rasterize_poly_refinement
    # - lowest refinement_level precedes
    expected = np.array(
        [
            [3, 3, 3, 3, 3, 3],
            [3, 2, 2, 2, 3, 3],
            [3, 2, 1, 1, 1, 3],
            [3, 2, 1, 1, 1, 3],
            [3, 2, 1, 1, 1, 3],
            [3, 2, 1, 1, 1, 3],
            [3, 3, 1, 1, 1, 3],
            [3, 3, 3, 3, 3, 3],
        ]
    )
    assert_array_equal(actual, expected)
