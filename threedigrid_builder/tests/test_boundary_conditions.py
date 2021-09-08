from numpy.testing import assert_array_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import BoundaryConditions1D
from threedigrid_builder.grid import Grid

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


def test_apply(grid1d):
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


def test_connection_node_does_not_exist(grid1d):
    boundary_conditions_1d = BoundaryConditions1D(id=[1], connection_node_id=[6])
    with pytest.raises(
        SchematisationError, match=r".*\[1\] refer to missing connection nodes \[6\].*"
    ):
        boundary_conditions_1d.apply(grid1d)


def test_connection_node_has_no_lines(grid1d):
    boundary_conditions_1d = BoundaryConditions1D(id=[1], connection_node_id=[5])
    grid1d.lines = grid1d.lines[[1, 2]]
    with pytest.raises(
        SchematisationError,
        match=r".*\[1\] refer to connection nodes that have no attached objects \(\[5\]\).*",
    ):
        boundary_conditions_1d.apply(grid1d)


def test_connection_node_has_too_many_lines(grid1d):
    boundary_conditions_1d = BoundaryConditions1D(id=[1], connection_node_id=[1])
    with pytest.raises(
        SchematisationError,
        match=r".*one attached object is allowed for connection nodes \[1\].*",
    ):
        boundary_conditions_1d.apply(grid1d)


def test_boundary_conditions_too_close(grid1d):
    boundary_conditions_1d = BoundaryConditions1D(id=[1, 2], connection_node_id=[1, 5])
    grid1d.lines = grid1d.lines[[0]]
    with pytest.raises(SchematisationError, match=r".*\[1, 2\] are too close to.*"):
        boundary_conditions_1d.apply(grid1d)
