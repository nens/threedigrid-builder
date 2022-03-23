from numpy.testing import assert_array_equal
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import GridRefinements
from threedigrid_builder.grid import QuadTree

import itertools
import numpy as np
import pygeos
import pytest


@pytest.fixture
def subgrid_meta():
    width = 20
    height = 16
    mask = np.ones((width, height), dtype=np.int16, order="F")
    mask[15:, :] = 0
    return {
        "pixel_size": 0.5,
        "width": width,
        "height": height,
        "bbox": (10.0, 10.0, 20.0, 18.0),
        "area_mask": mask,
    }


@pytest.fixture
def quadtree_no_refinement(subgrid_meta):
    return QuadTree(
        subgrid_meta=subgrid_meta,
        num_refine_levels=3,
        min_gridsize=1.0,
        use_2d_flow=True,
        refinements=None,
    )


@pytest.fixture
def quadtree_no_2d_flow(subgrid_meta):
    return QuadTree(
        subgrid_meta=subgrid_meta,
        num_refine_levels=3,
        min_gridsize=1.0,
        use_2d_flow=False,
        refinements=None,
    )


@pytest.fixture
def quadtree_line_refinement(subgrid_meta):
    refinement = GridRefinements(
        id=np.array([1]),
        refinement_level=np.array([1]),
        the_geom=pygeos.linestrings([[[15.0, 17.5], [17.5, 17.5], [17.5, 14.5]]]),
    )
    return QuadTree(
        subgrid_meta=subgrid_meta,
        num_refine_levels=3,
        min_gridsize=1.0,
        use_2d_flow=True,
        refinements=refinement,
    )


@pytest.fixture
def quadtree_poly_refinement(subgrid_meta):
    refinement = GridRefinements(
        id=np.array([1]),
        refinement_level=np.array([2]),
        the_geom=np.array([pygeos.box(12.4, 11.4, 13.9, 13.4)]),
    )
    return QuadTree(
        subgrid_meta=subgrid_meta,
        num_refine_levels=3,
        min_gridsize=1.0,
        use_2d_flow=True,
        refinements=refinement,
    )


def test_quadtree_no_refinement(quadtree_no_refinement):
    assert quadtree_no_refinement.kmax == 3
    assert np.size(quadtree_no_refinement.mmax) == 3
    assert quadtree_no_refinement.mmax[2] == 3
    assert quadtree_no_refinement.nmax[2] == 2
    assert quadtree_no_refinement.dx[0] == 1.0
    expected_lg = np.array(
        [
            [3, 3, 3, 3, 3, 3, 3, 3, -99, -99, -99, -99],
            [3, 3, 3, 3, 3, 3, 3, 3, -99, -99, -99, -99],
            [3, 3, 3, 3, 3, 3, 3, 3, -99, -99, -99, -99],
            [3, 3, 3, 3, 3, 3, 3, 3, -99, -99, -99, -99],
            [3, 3, 3, 3, 3, 3, 3, 3, -99, -99, -99, -99],
            [3, 3, 3, 3, 3, 3, 3, 3, -99, -99, -99, -99],
            [3, 3, 3, 3, 3, 3, 3, 3, -99, -99, -99, -99],
            [3, 3, 3, 3, 3, 3, 3, 3, -99, -99, -99, -99],
        ],
        dtype=np.int32,
    )
    assert_array_equal(quadtree_no_refinement.lg, expected_lg[::-1].T)


def test_quadtree_no_even_pixels(subgrid_meta):
    with pytest.raises(SchematisationError, match=r".*not contain an even number.*"):
        QuadTree(
            subgrid_meta=subgrid_meta,
            num_refine_levels=1,
            min_gridsize=1.5,
            use_2d_flow=True,
            refinements=None,
        )


def test_nodes_from_quadtree_line(quadtree_line_refinement, subgrid_meta):
    nodes, lines = quadtree_line_refinement.get_nodes_lines(
        subgrid_meta["area_mask"], itertools.count(start=0), itertools.count(start=0)
    )
    coordinates = np.array(
        [
            [12.0, 12.0],
            [11.0, 15.0],
            [11.0, 17.0],
            [13.0, 15.0],
            [13.0, 17.0],
            [15.0, 11.0],
            [15.0, 13.0],
            [15.0, 15.0],
            [17.0, 11.0],
            [17.0, 13.0],
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
            [17.5, 17.5],
        ]
    )

    pixel_coords_check = np.array(
        [[0, 0, 8, 8], [0, 12, 4, 16], [12, 0, 16, 4], [14, 14, 16, 16]]
    )
    assert_array_equal(nodes.coordinates, coordinates)
    assert_array_equal(nodes.pixel_coords[(0, 2, 8, 21), :], pixel_coords_check)


def test_lines_from_quadtree_line(quadtree_line_refinement, subgrid_meta):
    nodes, lines = quadtree_line_refinement.get_nodes_lines(
        subgrid_meta["area_mask"], itertools.count(start=0), itertools.count(start=0)
    )
    line = np.array(
        [
            [0, 5],
            [0, 6],
            [1, 3],
            [2, 4],
            [3, 7],
            [4, 10],
            [4, 11],
            [5, 8],
            [6, 9],
            [7, 14],
            [7, 15],
            [10, 12],
            [11, 13],
            [12, 16],
            [13, 17],
            [14, 18],
            [15, 19],
            [16, 20],
            [17, 21],
            [0, 1],
            [0, 3],
            [1, 2],
            [3, 4],
            [5, 6],
            [6, 7],
            [7, 10],
            [7, 12],
            [8, 9],
            [9, 14],
            [9, 18],
            [10, 11],
            [12, 13],
            [14, 15],
            [15, 16],
            [16, 17],
            [18, 19],
            [19, 20],
            [20, 21],
        ]
    )
    # fmt: off
    lik = np.array(
        [2, 2, 2, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
         2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], dtype=np.int32
    )
    lim = np.array(
        [3, 3, 1, 1, 2, 5, 5, 3, 3, 7, 7, 5, 5, 6, 6, 7, 7, 7, 7, 1, 2, 1, 2, 3,
         3, 5, 6, 4, 7, 8, 5, 6, 7, 7, 7, 8, 8, 8], dtype=np.int32
    )
    lin = np.array(
        [1, 2, 3, 4, 3, 7, 8, 1, 2, 5, 6, 7, 8, 7, 8, 5, 6, 7, 8, 3, 3, 3, 3, 1,
         2, 7, 7, 1, 5, 5, 7, 7, 5, 6, 7, 5, 6, 7], dtype=np.int32
    )
    # fmt: on

    cross_pix_coords = np.array(
        [
            [8, 0, 8, 4],
            [8, 4, 8, 8],
            [4, 8, 4, 12],
            [8, 12, 10, 12],
            [14, 8, 16, 8],
            [14, 10, 16, 10],
        ],
        dtype=np.int32,
    )

    assert_array_equal(lines.line, line)
    assert_array_equal(lines.lik, lik)
    assert_array_equal(lines.lim, lim)
    assert_array_equal(lines.lin, lin)
    assert_array_equal(
        lines.cross_pix_coords[(0, 1, 2, 25, 29, 35), :], cross_pix_coords
    )


def test_nodes_from_quadtree_poly(quadtree_poly_refinement, subgrid_meta):
    nodes, lines = quadtree_poly_refinement.get_nodes_lines(
        subgrid_meta["area_mask"], itertools.count(start=0), itertools.count(start=0)
    )
    expected_coordinates = [
        [12.0, 16.0],
        [16.0, 12.0],
        [16.0, 16.0],
        [11.0, 11.0],
        [11.0, 13.0],
        [13.0, 11.0],
        [13.0, 13.0],
    ]

    expected_pixel_coords = [
        [0, 8, 8, 16],
        [8, 0, 16, 8],
        [8, 8, 16, 16],
        [0, 0, 4, 4],
        [0, 4, 4, 8],
        [4, 0, 8, 4],
        [4, 4, 8, 8],
    ]
    assert_array_equal(nodes.coordinates, expected_coordinates)
    assert_array_equal(nodes.pixel_coords, expected_pixel_coords)


def test_lines_from_quadtree_poly(quadtree_poly_refinement, subgrid_meta):
    nodes, lines = quadtree_poly_refinement.get_nodes_lines(
        subgrid_meta["area_mask"], itertools.count(start=0), itertools.count(start=0)
    )
    line = [
        [0, 2],
        [5, 1],
        [6, 1],
        [3, 5],
        [4, 6],
        [4, 0],
        [6, 0],
        [1, 2],
        [3, 4],
        [5, 6],
    ]
    lik = [3, 2, 2, 2, 2, 2, 2, 3, 2, 2]
    lim = [1, 2, 2, 1, 1, 1, 2, 2, 1, 2]
    lin = [2, 1, 2, 1, 2, 2, 2, 1, 1, 1]

    cross_pix_coords = [
        [8, 8, 8, 16],
        [8, 0, 8, 4],
        [8, 4, 8, 8],
        [4, 0, 4, 4],
        [4, 4, 4, 8],
        [0, 8, 4, 8],
        [4, 8, 8, 8],
        [8, 8, 16, 8],
        [0, 4, 4, 4],
        [4, 4, 8, 4],
    ]
    assert_array_equal(lines.line, line)
    assert_array_equal(lines.lik, lik)
    assert_array_equal(lines.lim, lim)
    assert_array_equal(lines.lin, lin)
    assert_array_equal(lines.cross_pix_coords, cross_pix_coords)


def test_no_2d_flow(quadtree_no_2d_flow, subgrid_meta):
    nodes, lines = quadtree_no_2d_flow.get_nodes_lines(
        subgrid_meta["area_mask"],
        itertools.count(start=0),
        itertools.count(start=0),
    )

    assert len(lines) == 0
    assert len(nodes) == 4
