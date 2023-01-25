import itertools
from unittest import mock

import numpy as np
import pytest
import shapely
from numpy.testing import assert_array_equal

from threedigrid_builder.base import Lines, Nodes
from threedigrid_builder.constants import (
    BoundaryType,
    CalculationType,
    ContentType,
    LineType,
    NodeType,
)
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import BoundaryConditions1D, BoundaryConditions2D, Grid


@pytest.fixture
def grid1d():
    return Grid(
        nodes=Nodes(
            id=[1, 2, 3, 4],
            calculation_type=CalculationType.ISOLATED,
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            node_type=NodeType.NODE_1D_NO_STORAGE,
            content_pk=[1, 5, 9, 10],
        ),
        lines=Lines(
            id=[1, 2, 3],
            line=[(2, 1), (1, 3), (1, 4)],
            kcu=LineType.LINE_1D_ISOLATED,
        ),
    )


def test_1d_apply(grid1d):
    boundary_conditions_1d = BoundaryConditions1D(
        id=[10, 22],
        boundary_type=[2, 3],
        connection_node_id=[9, 5],
    )
    boundary_conditions_1d.apply(grid1d)

    # note that apply does not reorder the nodes, this happens later at grid.sort()
    assert_array_equal(grid1d.nodes.boundary_id, [-9999, 22, 10, -9999])
    assert_array_equal(grid1d.nodes.boundary_type, [-9999, 3, 2, -9999])
    assert_array_equal(grid1d.nodes.node_type[[0, 3]], NodeType.NODE_1D_NO_STORAGE)
    assert_array_equal(grid1d.nodes.node_type[[1, 2]], NodeType.NODE_1D_BOUNDARIES)
    assert_array_equal(grid1d.nodes.calculation_type[[0, 3]], CalculationType.ISOLATED)
    assert_array_equal(
        grid1d.nodes.calculation_type[[1, 2]], CalculationType.BOUNDARY_NODE
    )
    assert_array_equal(grid1d.lines.is_1d_boundary[[0, 1]], 1)
    assert_array_equal(grid1d.lines.is_1d_boundary[[2]], -9999)


def test_1d_connection_node_does_not_exist(grid1d):
    boundary_conditions_1d = BoundaryConditions1D(id=[1], connection_node_id=[6])
    with pytest.raises(
        SchematisationError, match=r".*\[1\] refer to missing connection nodes \[6\].*"
    ):
        boundary_conditions_1d.apply(grid1d)


def test_1d_connection_node_has_no_lines(grid1d):
    boundary_conditions_1d = BoundaryConditions1D(id=[1], connection_node_id=[5])
    grid1d.lines = grid1d.lines[[1, 2]]
    with pytest.raises(
        SchematisationError,
        match=r".*\[1\] refer to connection nodes that have no attached objects \(\[5\]\).*",
    ):
        boundary_conditions_1d.apply(grid1d)


def test_1d_connection_node_has_too_many_lines(grid1d):
    boundary_conditions_1d = BoundaryConditions1D(id=[1], connection_node_id=[1])
    with pytest.raises(
        SchematisationError,
        match=r".*one attached object is allowed for connection nodes \[1\].*",
    ):
        boundary_conditions_1d.apply(grid1d)


def test_1d_boundary_conditions_too_close(grid1d):
    boundary_conditions_1d = BoundaryConditions1D(id=[1, 2], connection_node_id=[1, 5])
    grid1d.lines = grid1d.lines[[0]]
    with pytest.raises(SchematisationError, match=r".*\[1, 2\] are too close to.*"):
        boundary_conditions_1d.apply(grid1d)


@pytest.fixture
def grid2d():
    # A grid as follows, without lines.
    #
    #       + - + - - - +
    #       | 5 |       |
    #   + - + - +   2   |
    #   | 3 | 4 |       |
    #   + - + - + - + - +
    #   |       | 7 | 8 |
    #   |   1   + - + - +
    #   |       | 6 |
    #   + - - - + - +
    #
    pixel_size = 5.0
    min_cell_size = 20.0  # the size of the smallest cell
    origin = 46.0, 47.0  # coordinate of the bottomleft corner
    nodk = np.array([2, 2, 1, 1, 1, 1, 1, 1])
    nodm = np.array([1, 2, 1, 2, 2, 3, 3, 4])
    nodn = np.array([1, 2, 3, 3, 4, 1, 2, 2])
    _size = 2 ** (nodk - 1) * min_cell_size
    bounds = np.empty((len(nodk), 4))
    bounds[:, 0] = (nodm - 1) * _size + origin[0]
    bounds[:, 1] = (nodn - 1) * _size + origin[1]
    bounds[:, 2] = bounds[:, 0] + _size
    bounds[:, 3] = bounds[:, 1] + _size
    pixel_coords = np.empty((len(nodk), 4), dtype=int)
    pixel_coords[:, 0] = (nodm - 1) * _size / pixel_size
    pixel_coords[:, 1] = (nodn - 1) * _size / pixel_size
    pixel_coords[:, 2] = pixel_coords[:, 0] + _size / pixel_size
    pixel_coords[:, 3] = pixel_coords[:, 1] + _size / pixel_size
    nodes = Nodes(
        id=range(1, 9),
        node_type=NodeType.NODE_2D_OPEN_WATER,
        bounds=bounds,
        pixel_coords=pixel_coords,
        nodk=nodk,
        nodm=nodm,
        nodn=nodn,
    )
    grid = Grid(nodes=nodes, lines=Lines(id=[]))
    grid.quadtree = mock.Mock()
    grid.quadtree.quad_idx = np.array(
        [
            [0, 5, 2, 2],
            [3, 4, 2, 2],
            [1, 1, 7, 8],
            [1, 1, 6, 0],
        ]
    )[::-1].T
    return grid


# some expected boundary cells
bottom_1 = [46, 7, 86, 47]
bottom_6 = [86, 27, 106, 47]
left_1 = [6, 47, 46, 87]
left_3 = [26, 87, 46, 107]
top_2 = [86, 127, 126, 167]
top_5 = [66, 127, 86, 147]
right_2 = [126, 87, 166, 127]
right_8 = [126, 67, 146, 87]

SOUTH = LineType.LINE_2D_BOUNDARY_SOUTH
NORTH = LineType.LINE_2D_BOUNDARY_NORTH
EAST = LineType.LINE_2D_BOUNDARY_EAST
WEST = LineType.LINE_2D_BOUNDARY_WEST


@pytest.mark.parametrize(
    "bc_coords, bounds, boundary_id, nodm, nodn, kcu, line, cross_pix_coords",
    [
        (  # Boundary condition exactly at the bottom
            [[(46, 47), (106, 47)]],
            (bottom_1, bottom_6),
            0,
            [1, 3],
            [0, 0],
            SOUTH,
            [(9, 1), (10, 6)],
            [(0, 0, 8, 0), (8, 0, 12, 0)],
        ),
        (  # Boundary condition exactly at the left
            [[(46, 47), (46, 107)]],
            (left_1, left_3),
            0,
            [0, 0],
            [1, 3],
            WEST,
            [(9, 1), (10, 3)],
            [(0, 0, 0, 8), (0, 8, 0, 12)],
        ),
        (  # Boundary condition exactly at the top
            [[(66, 127), (126, 127)]],
            (top_2, top_5),
            0,
            [2, 2],
            [3, 5],
            NORTH,
            [(2, 9), (5, 10)],
            [(8, 16, 16, 16), (4, 16, 8, 16)],
        ),
        (  # Boundary condition exactly at the right
            [[(126, 67), (126, 127)]],
            (right_2, right_8),
            0,
            [3, 5],
            [2, 2],
            EAST,
            [(2, 9), (8, 10)],
            [(16, 8, 16, 16), (16, 4, 16, 8)],
        ),
        (  # Boundary condition exactly at the bottom and then one left (with overlap)
            [[(46, 47), (106, 47)], [(46, 47), (46, 107)]],
            (bottom_1, bottom_6, left_3),
            (0, 0, 1),
            [1, 3, 0],
            [0, 0, 3],
            [SOUTH, SOUTH, WEST],
            [(9, 1), (10, 6), (11, 3)],
            [(0, 0, 8, 0), (8, 0, 12, 0), (0, 8, 0, 12)],
        ),
        (  # Boundary condition exactly at the left and then one bottom (with overlap)
            [[(46, 47), (46, 107)], [(46, 47), (106, 47)]],
            (left_1, left_3, bottom_6),
            (0, 0, 1),
            [0, 0, 3],
            [1, 3, 0],
            [WEST, WEST, SOUTH],
            [(9, 1), (10, 3), (11, 6)],
            [(0, 0, 0, 8), (0, 8, 0, 12), (8, 0, 12, 0)],
        ),
        (  # Boundary condition just outside of the model, bottom
            [[(46, 40), (106, 41)]],
            (bottom_1, bottom_6),
            0,
            [1, 3],
            [0, 0],
            SOUTH,
            [(9, 1), (10, 6)],
            [(0, 0, 8, 0), (8, 0, 12, 0)],
        ),
        (  # Boundary condition just outside of the model, left
            [[(40, 47), (41, 107)]],
            (left_1, left_3),
            0,
            [0, 0],
            [1, 3],
            WEST,
            [(9, 1), (10, 3)],
            [(0, 0, 0, 8), (0, 8, 0, 12)],
        ),
        (  # Boundary condition just outside of the model, top
            [[(66, 130), (126, 132)]],
            (top_2, top_5),
            0,
            [2, 2],
            [3, 5],
            NORTH,
            [(2, 9), (5, 10)],
            [(8, 16, 16, 16), (4, 16, 8, 16)],
        ),
        (  # Boundary condition just outside of the model, right
            [[(127, 67), (126, 127)]],
            (right_2, right_8),
            0,
            [3, 5],
            [2, 2],
            EAST,
            [(2, 9), (8, 10)],
            [(16, 8, 16, 16), (16, 4, 16, 8)],
        ),
    ],
)
def test_2d_boundary_condition(
    grid2d,
    bc_coords,
    bounds,
    nodm,
    nodn,
    kcu,
    boundary_id,
    line,
    cross_pix_coords,
):
    boundary_conditions_2d = BoundaryConditions2D(
        id=range(len(bc_coords)),
        boundary_type=[3] * len(bc_coords),
        the_geom=shapely.linestrings(bc_coords),
    )

    nodes, lines = boundary_conditions_2d.get_nodes_and_lines(
        grid2d.nodes,
        grid2d.cell_tree,
        grid2d.quadtree,
        itertools.count(9),
        itertools.count(),
    )

    assert_array_equal(nodes.node_type, NodeType.NODE_2D_BOUNDARIES)
    assert_array_equal(nodes.boundary_id, boundary_id)
    assert_array_equal(nodes.boundary_type, 3)
    assert_array_equal(nodes.bounds, bounds)
    assert_array_equal(nodes.pixel_coords, -9999)
    assert_array_equal(nodes.nodm, nodm)
    assert_array_equal(nodes.nodn, nodn)
    assert_array_equal(lines.kcu, kcu)
    assert_array_equal(lines.line, line)
    assert_array_equal(lines.cross_pix_coords, cross_pix_coords)


@pytest.mark.parametrize(
    "bc_coords, expected_message",
    [
        ([[(46, 0), (106, 0)]], r".*does not touch any edge cell."),
        (
            [[(46, 77), (106, 77)]],
            r".*different edge coordinates \(y=\[47.0, 67.0\]\).",
        ),
        (
            [[(76, 47), (76, 107)]],
            r".*different edge coordinates \(x=\[46.0, 66.0\]\).",
        ),
    ],
)
def test_2d_boundary_condition_err(grid2d, bc_coords, expected_message):
    boundary_conditions_2d = BoundaryConditions2D(
        id=range(len(bc_coords)),
        boundary_type=[3] * len(bc_coords),
        the_geom=shapely.linestrings(bc_coords),
    )

    with pytest.raises(SchematisationError, match=expected_message):
        boundary_conditions_2d.get_nodes_and_lines(
            grid2d.nodes,
            grid2d.cell_tree,
            grid2d.quadtree,
            itertools.count(),
            itertools.count(),
        )


@pytest.fixture
def grid2d_gw():
    nodes = Nodes(
        id=range(4),
        node_type=([NodeType.NODE_2D_OPEN_WATER] * 2)
        + ([NodeType.NODE_2D_GROUNDWATER] * 2),
        bounds=[(2.0, 3.0, 7.0, 8.0), (7.0, 3.0, 12.0, 8.0)] * 2,
        pixel_coords=[(0, 0, 10, 10), (10, 0, 20, 10)] * 2,
        nodk=[1, 1, 1, 1],
        nodm=[1, 2, 1, 2],
        nodn=[1, 1, 1, 1],
    )
    grid = Grid(nodes=nodes, lines=Lines(id=[]))
    grid.quadtree = mock.Mock()
    grid.quadtree.quad_idx = np.array([[1, 2]])[::-1].T
    return grid


@pytest.mark.parametrize(
    "boundary_type,node_type,kcu,line",
    [
        (
            BoundaryType.WATERLEVEL,
            NodeType.NODE_2D_BOUNDARIES,
            LineType.LINE_2D_BOUNDARY_EAST,
            (1, 9),
        ),
        (
            BoundaryType.DISCHARGE,
            NodeType.NODE_2D_BOUNDARIES,
            LineType.LINE_2D_BOUNDARY_EAST,
            (1, 9),
        ),
        (
            BoundaryType.GROUNDWATERLEVEL,
            NodeType.NODE_2D_GROUNDWATER_BOUNDARIES,
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_EAST,
            (3, 9),
        ),
        (
            BoundaryType.GROUNDWATERDISCHARGE,
            NodeType.NODE_2D_GROUNDWATER_BOUNDARIES,
            LineType.LINE_2D_GROUNDWATER_BOUNDARY_EAST,
            (3, 9),
        ),
    ],
)
def test_2d_boundary_condition_types(grid2d_gw, boundary_type, node_type, kcu, line):
    boundary_conditions_2d = BoundaryConditions2D(
        id=[1],
        boundary_type=[boundary_type],
        the_geom=[shapely.linestrings([(12.0, 3.0), (12.0, 8.0)])],
    )
    nodes, lines = boundary_conditions_2d.get_nodes_and_lines(
        grid2d_gw.nodes,
        grid2d_gw.cell_tree,
        grid2d_gw.quadtree,
        itertools.count(9),
        itertools.count(7),
    )

    assert_array_equal(nodes.id, [9])
    assert_array_equal(nodes.node_type, [node_type])
    assert_array_equal(nodes.boundary_id, [1])
    assert_array_equal(nodes.boundary_type, [boundary_type])
    assert_array_equal(nodes.bounds, [(12.0, 3.0, 17, 8.0)])
    assert_array_equal(nodes.pixel_coords, -9999)
    assert_array_equal(nodes.nodm, [3])
    assert_array_equal(nodes.nodn, [1])
    assert_array_equal(lines.id, [7])
    assert_array_equal(lines.kcu, [kcu])
    assert_array_equal(lines.line, [line])
    assert_array_equal(lines.cross_pix_coords, [(20, 0, 20, 10)])


def test_2d_boundary_condition_combined(grid2d_gw):
    boundary_conditions_2d = BoundaryConditions2D(
        id=[1, 2],
        boundary_type=[BoundaryType.WATERLEVEL, BoundaryType.GROUNDWATERLEVEL],
        the_geom=[shapely.linestrings([(12.0, 3.0), (12.0, 8.0)])] * 2,
    )

    nodes, lines = boundary_conditions_2d.get_nodes_and_lines(
        grid2d_gw.nodes,
        grid2d_gw.cell_tree,
        grid2d_gw.quadtree,
        itertools.count(9),
        itertools.count(7),
    )

    assert_array_equal(nodes.id, [9, 10])
    assert_array_equal(
        nodes.node_type,
        [NodeType.NODE_2D_BOUNDARIES, NodeType.NODE_2D_GROUNDWATER_BOUNDARIES],
    )
    assert_array_equal(nodes.boundary_id, [1, 2])
    assert_array_equal(
        nodes.boundary_type, [BoundaryType.WATERLEVEL, BoundaryType.GROUNDWATERLEVEL]
    )
    assert_array_equal(nodes.bounds, [(12.0, 3.0, 17, 8.0)] * 2)
    assert_array_equal(nodes.pixel_coords, -9999)
    assert_array_equal(nodes.nodm, 3)
    assert_array_equal(nodes.nodn, 1)
    assert_array_equal(lines.id, [7, 8])
    assert_array_equal(
        lines.kcu,
        [LineType.LINE_2D_BOUNDARY_EAST, LineType.LINE_2D_GROUNDWATER_BOUNDARY_EAST],
    )
    assert_array_equal(lines.line, [(1, 9), (3, 10)])
    assert_array_equal(lines.cross_pix_coords, [(20, 0, 20, 10)] * 2)
