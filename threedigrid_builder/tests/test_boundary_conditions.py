from numpy.testing import assert_array_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import BoundaryConditions1D
from threedigrid_builder.grid import BoundaryConditions2D
from threedigrid_builder.grid import Grid
from unittest import mock

import itertools
import numpy as np
import pygeos
import pytest


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


def test_1d_pply(grid1d):
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
    assert_array_equal(grid1d.lines.kcu[[0, 1]], LineType.LINE_1D_BOUNDARY)
    assert_array_equal(grid1d.lines.kcu[[2]], LineType.LINE_1D_ISOLATED)


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
    nodes = Nodes(
        id=range(1, 9),
        node_type=NodeType.NODE_2D_OPEN_WATER,
        bounds=bounds,
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


@pytest.mark.parametrize(
    "bc_coords, expected_bounds",
    [
        ([[(46, 47), (106, 47)]], (bottom_1, bottom_6)),  # exact
        ([[(50, 47), (120, 66)]], (bottom_1, bottom_6)),  # inexact
        ([[(46, 47), (46, 107)]], (left_1, left_3)),  # exact
        ([[(50, 50), (65, 113)]], (left_1, left_3)),  # inexact
        ([[(66, 127), (126, 127)]], (top_2, top_5)),  # exact
        ([[(56, 125), (100, 110)]], (top_2, top_5)),  # inexact
        ([[(126, 67), (126, 127)]], (right_2, right_8)),  # exact
        ([[(127, 70), (110, 140)]], (right_2, right_8)),  # inexact
        ([[(46, 47), (106, 47)], [(46, 47), (46, 107)]], (bottom_1, bottom_6, left_3)),
        ([[(46, 47), (46, 107)], [(46, 47), (106, 47)]], (left_1, left_3, bottom_6)),
    ],
)
def test_2d_boundary_condition(grid2d, bc_coords, expected_bounds):
    boundary_conditions_2d = BoundaryConditions2D(
        id=range(len(bc_coords)),
        boundary_type=[3] * len(bc_coords),
        the_geom=pygeos.linestrings(bc_coords),
    )

    nodes = boundary_conditions_2d.get_nodes(
        grid2d.nodes, grid2d.cell_tree, grid2d.quadtree, itertools.count()
    )
    assert_array_equal(nodes.bounds, expected_bounds)


@pytest.mark.parametrize(
    "bc_coords, expected_message",
    [
        ([[(46, 30), (106, 30)]], r".*does not touch any edge cell."),
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
        the_geom=pygeos.linestrings(bc_coords),
    )

    with pytest.raises(SchematisationError, match=expected_message):
        boundary_conditions_2d.get_nodes(
            grid2d.nodes, grid2d.cell_tree, grid2d.quadtree, itertools.count()
        )
