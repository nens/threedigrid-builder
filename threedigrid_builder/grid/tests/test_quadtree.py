from numpy.testing import assert_array_equal
from threedigrid_builder.base import Nodes
from threedigrid_builder.grid import QuadTree
from threedigrid_builder.grid import GridRefinements

import numpy as np
import pygeos
import pytest


@pytest.fixture
def subgrid_meta():
    width = 20
    height = 16
    mask = np.ones((width, height), dtype=np.int32, order="F")
    mask[15:, :] = 0
    return{
        "pixel_size": 0.5,
        "width": width,
        "height": height,
        "bbox": (10., 10., 20., 18.),
        "area_mask": mask,
    }


@pytest.fixture
def quadtree_no_refinement(subgrid_meta):
    return QuadTree(
        subgrid_meta=subgrid_meta,
        num_refine_levels=3,
        min_gridsize=1.0,
        refinements=None,
    )


@pytest.fixture
def quadtree_line_refinement(subgrid_meta):
    refinement = GridRefinements(
        id=np.array([1]),
        refinement_level=np.array([1]),
        the_geom=np.array(
            [pygeos.linestrings([[15., 17.5], [17.5, 17.5], [17.5, 14.5]])]
        )
    )
    return QuadTree(
        subgrid_meta=subgrid_meta,
        num_refine_levels=3,
        min_gridsize=1.0,
        refinements=refinement,
    )


@pytest.fixture
def quadtree_poly_refinement(subgrid_meta):
    refinement = GridRefinements(
        id=np.array([1]),
        refinement_level=np.array([2]),
        the_geom=np.array(
            [pygeos.box(15., 17.5, 17.5, 14.5)]
        )
    )
    return QuadTree(
        subgrid_meta=subgrid_meta,
        num_refine_levels=3,
        min_gridsize=1.0,
        refinements=refinement,
    )


def test_quadtree_no_refinement(quadtree_no_refinement):
    assert quadtree_no_refinement.kmax == 3
    assert np.size(quadtree_no_refinement.mmax) == 3
    assert quadtree_no_refinement.mmax[2] == 3
    assert quadtree_no_refinement.nmax[2] == 2
    assert quadtree_no_refinement.dx[0] == 1.0
    assert_array_equal(
        quadtree_no_refinement.lg,
        np.array(
            [[3, 3, 3, 3, 3, 3, 3, 3],
             [3, 3, 3, 3, 3, 3, 3, 3],
             [3, 3, 3, 3, 3, 3, 3, 3],
             [3, 3, 3, 3, 3, 3, 3, 3],
             [3, 3, 3, 3, 3, 3, 3, 3],
             [3, 3, 3, 3, 3, 3, 3, 3],
             [3, 3, 3, 3, 3, 3, 3, 3],
             [3, 3, 3, 3, 3, 3, 3, 3],
             [-99, -99, -99, -99, -99, -99, -99, -99],
             [-99, -99, -99, -99, -99, -99, -99, -99],
             [-99, -99, -99, -99, -99, -99, -99, -99],
             [-99, -99, -99, -99, -99, -99, -99, -99]],
            dtype=np.int32, order='F'
        )
    )


def test_quadtree_line_refinement(quadtree_line_refinement):

    assert quadtree_line_refinement.kmax == 3
    assert np.size(quadtree_line_refinement.mmax) == 3
    assert quadtree_line_refinement.mmax[2] == 3
    assert quadtree_line_refinement.nmax[2] == 2
    assert quadtree_line_refinement.dx[0] == 1.0
    assert_array_equal(
        quadtree_line_refinement.lg,
        np.array(
            [[3, 3, 3, 3, 2, 2, 2, 2],
             [3, 3, 3, 3, 2, 2, 2, 2],
             [3, 3, 3, 3, 2, 2, 2, 2],
             [3, 3, 3, 3, 2, 2, 2, 2],
             [2, 2, 2, 2, 2, 2, 1, 1],
             [2, 2, 2, 2, 2, 2, 1, 1],
             [2, 2, 2, 2, 1, 1, 1, 1],
             [2, 2, 2, 2, 1, 1, 1, 1],
             [-99, -99, -99, -99, -99, -99, -99, -99],
             [-99, -99, -99, -99, -99, -99, -99, -99],
             [-99, -99, -99, -99, -99, -99, -99, -99],
             [-99, -99, -99, -99, -99, -99, -99, -99]],
            dtype=np.int32, order='F'
        )
    )


def test_quadtree_poly_refinement(quadtree_poly_refinement):

    assert quadtree_poly_refinement.kmax == 3
    assert np.size(quadtree_poly_refinement.mmax) == 3
    assert quadtree_poly_refinement.mmax[2] == 3
    assert quadtree_poly_refinement.nmax[2] == 2
    assert quadtree_poly_refinement.dx[0] == 1.0
    assert_array_equal(
        quadtree_poly_refinement.lg,
        np.array(
            [[3, 3, 3, 3, 3, 3, 3, 3],
             [3, 3, 3, 3, 3, 3, 3, 3],
             [3, 3, 3, 3, 3, 3, 3, 3],
             [3, 3, 3, 3, 3, 3, 3, 3],
             [3, 3, 3, 3, 2, 2, 2, 2],
             [3, 3, 3, 3, 2, 2, 2, 2],
             [3, 3, 3, 3, 2, 2, 2, 2],
             [3, 3, 3, 3, 2, 2, 2, 2],
             [-99, -99, -99, -99, -99, -99, -99, -99],
             [-99, -99, -99, -99, -99, -99, -99, -99],
             [-99, -99, -99, -99, -99, -99, -99, -99],
             [-99, -99, -99, -99, -99, -99, -99, -99]],
            dtype=np.int32, order='F'
        )
    )


def test_nodes_from_quadtree(quadtree_line_refinement, subgrid_meta):
    nodes, lines = quadtree_line_refinement.get_nodes_lines(
        subgrid_meta["area_mask"]
    )
    coordinates = np.array(
        [[12., 12.],
         [11., 15.],
         [11., 17.],
         [13., 15.],
         [13., 17.],
         [15., 11.],
         [15., 13.],
         [15., 15.],
         [17., 11.],
         [17., 13.],
         [14.5, 16.5],
         [14.5, 17.5],
         [15.5, 16.5],
         [15.5, 17.5],
         [16.5, 14.5],
         [16.5, 15.5],
         [16.5, 16.5],
         [16.5, 17.5],
         [17.5, 14.5],
         [17.5, 15.5],
         [17.5, 16.5],
         [17.5, 17.5]]
    )
    assert_array_equal(nodes.coordinates, coordinates)


def test_lines_from_quadtree(quadtree_line_refinement, subgrid_meta):

    nodes, lines = quadtree_line_refinement.get_nodes_lines(
        subgrid_meta["area_mask"]
    )

    line = np.array(
        [[0, 5], [0, 6], [1, 3], [2, 4], [3, 7], [4, 10], [4, 11], [5, 8], [6, 9],
         [7, 14], [7, 15], [10, 12], [11, 13], [12, 16], [13, 17], [14, 18],
         [15, 19], [16, 20], [17, 21],
         [0, 1], [0, 3], [1, 2], [3, 4], [5, 6], [6, 7], [7, 10], [7, 12],
         [8, 9], [9, 14], [9, 18], [10, 11], [12, 13], [14, 15], [15, 16],
         [16, 17], [18, 19], [19, 20], [20, 21]]
    )
    assert_array_equal(lines.line, line)
