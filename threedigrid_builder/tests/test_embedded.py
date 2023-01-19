from itertools import count

import pytest
import shapely
from numpy.testing import assert_almost_equal, assert_array_equal

from threedigrid_builder.base import Array, Lines, Nodes
from threedigrid_builder.constants import CalculationType, ContentType, NodeType
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import (
    ConnectionNodes,
    embed_linear_objects,
    embed_nodes,
    Grid,
)
from threedigrid_builder.grid.linear import BaseLinear

NODE_2D = NodeType.NODE_2D_OPEN_WATER
NODE_1D = NodeType.NODE_1D_STORAGE
EMBEDDED = CalculationType.EMBEDDED
ISOLATED = CalculationType.ISOLATED


class LinearObject:
    id: int
    the_geom: shapely.Geometry
    calculation_type: CalculationType
    connection_node_start_id: int
    connection_node_end_id: int
    display_name: str


class LinearObjects(Array[LinearObject], BaseLinear):
    content_type = ContentType.TYPE_V2_WINDSHIELD  # just pick one for the test


@pytest.fixture
def grid2d():
    return Grid(
        nodes=Nodes(
            id=[0, 1, 2, 3],
            node_type=NODE_2D,
            bounds=[(0, 0, 10, 10), (10, 0, 20, 10), (0, 10, 10, 20), (10, 10, 20, 20)],
        ),
        lines=Lines(id=[0], line=[[0, 1]]),
    )


def test_embed_nodes(grid2d):
    grid = grid2d + Grid(
        nodes=Nodes(
            id=[4, 5, 6, 7],
            node_type=NODE_1D,
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            calculation_type=[EMBEDDED, ISOLATED, EMBEDDED, EMBEDDED],
            coordinates=[(10, 2), (9, 9), (13, 3), (16, 7)],
        ),
        lines=Lines(id=[1, 2, 3], line=[(4, 5), (4, 6), (5, 6)]),
    )
    embedded = embed_nodes(grid, count(10))

    assert_array_equal(grid.nodes.id, [0, 1, 2, 3, 4])
    assert_array_equal(grid.nodes.node_type, [NODE_2D] * 4 + [NODE_1D])
    assert_array_equal(grid.lines.line, [(0, 1), (0, 4), (0, 1), (4, 1)])

    assert_array_equal(embedded.id, [10, 11, 12])  # new ids
    assert_array_equal(embedded.calculation_type, EMBEDDED)
    assert_array_equal(embedded.embedded_in, [0, 1, 1])


def test_embed_node_outside_2D(grid2d):
    grid = grid2d + Grid(
        nodes=Nodes(
            id=[4, 5],
            node_type=NODE_1D,
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            content_pk=[15, 16],
            calculation_type=EMBEDDED,
            coordinates=[(5, 22), (5, 5)],
        ),
        lines=Lines(id=[]),
    )

    with pytest.raises(SchematisationError, match=r".*\[15\] are outside the 2D cells"):
        embed_nodes(grid, count(10))


def test_embed_node_two_interconnected(grid2d):
    grid = grid2d + Grid(
        nodes=Nodes(
            id=[4, 5],
            node_type=NODE_1D,
            content_type=ContentType.TYPE_V2_CONNECTION_NODES,
            content_pk=[4, 5],
            calculation_type=EMBEDDED,
            coordinates=[(8, 8), (5, 5)],
        ),
        lines=Lines(id=[2], line=[(4, 5)]),
    )

    with pytest.raises(SchematisationError, match=r".*\[4, 5\] connect to.*"):
        embed_nodes(grid, count(10))


@pytest.fixture
def connection_nodes():
    # Used to map connection_node_start/end_id to an index (sequence id)
    return ConnectionNodes(
        id=[21, 25, 33, 42],
        the_geom=shapely.points([(0, 21), (1, 25), (2, 33), (3, 42)]),
    )


def test_embed_linear_objects_multiple(grid2d):
    linear_objects = LinearObjects(
        id=[0, 1, 2, 3],
        calculation_type=EMBEDDED,
        the_geom=[
            shapely.linestrings([(1, 1), (15, 1)]),
            shapely.linestrings([(11, 1), (15, 1)]),
            shapely.linestrings([(6, 7), (10, 10), (15, 10), (15, 15), (10, 15)]),
            shapely.linestrings([(7, 2), (18, 2), (18, 9), (18, 18)]),
        ],
    )

    nodes, lines_s1d, lines_ds1d_half = embed_linear_objects(
        linear_objects,
        grid2d.cell_tree,
        embedded_cutoff_threshold=1,
        embedded_node_id_counter=count(2),
    )

    assert_array_equal(nodes.id, [2])
    assert_array_equal(nodes.content_type, ContentType.TYPE_V2_WINDSHIELD)
    assert_array_equal(nodes.content_pk, [3])
    assert_array_equal(nodes.coordinates, [(18, 2)])
    assert_array_equal(nodes.calculation_type, EMBEDDED)
    assert_array_equal(nodes.s1d, [11])  # halfway the 2 velocity points (see below)
    assert_array_equal(nodes.embedded_in, [1])
    assert_almost_equal(lines_s1d, [9, 2, 5, 3, 19])
    assert_almost_equal(lines_ds1d_half, [9, 2, 5, 3, 8])


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize(
    "geometry,lines_s1d,embedded_in",
    [
        ([(5, 5), (9, 5)], [2], []),  # no cell crossing --> no nodes and 1 line
        ([(5, 5), (18, 5)], [5], []),  # horizontal, 1 crossing
        ([(5, 5), (5, 18)], [5], []),  # vertical, 1 crossing
        ([(5, 5), (5, 9), (15, 9)], [9], []),  # 1 crossing, more coords
        ([(5, 5), (18, 5), (18, 9)], [5], []),  # 1 crossing, more coords
        ([(5, 15), (5, 1), (18, 1)], [5, 19], [0]),  # corner, bottomleft
        ([(5, 5), (18, 5), (18, 13)], [5, 18], [1]),  # corner, bottomright
        ([(5, 5), (5, 14), (18, 14)], [5, 14], [2]),  # corner, bottomright
        ([(18, 5), (18, 13), (9, 13)], [5, 16], [3]),  # corner, topright
        ([(5, 15), (5, 2), (19, 2), (19, 15)], [5, 18, 35], [0, 1]),  # U
        ([(5, 2), (19, 2), (19, 15), (5, 15)], [5, 22, 36], [1, 3]),  # U, on left side
        ([(19, 2), (19, 15), (5, 15), (5, 2)], [8, 22, 32], [3, 2]),  # U, upside down
        ([(19, 15), (5, 15), (5, 2), (19, 2)], [9, 19, 32], [2, 0]),  # U, on right side
        ([(8, 1), (12, 4), (8, 7)], [2.5, 7.5], [1]),  # 0 - 1 - 0
        ([(8, 1), (12, 4), (8, 7), (12, 7)], [2.5, 7.5, 12.0], [1, 0]),  # 0 - 1 - 0 - 1
        ([(5, 5), (18, 5), (18, 10)], [5], []),  # 1 crossing, end at edge
        ([(5, 15), (5, 1), (18, 1), (18, 10)], [5, 19], [0]),  # corner, end at edge
        ([(6, 1), (10, 4), (6, 7)], [5.0], []),  # 0 - touch 1
        ([(14, 1), (10, 4), (14, 7)], [5.0], []),  # 1 - touch 0
        ([(12, 7), (6, 7), (10, 4), (6, 1)], [2], []),  # 2 - 0 - touch 1
        ([(6, 1), (10, 4), (6, 7), (12, 7)], [14], []),  # 0 - touch 1 - 2
        ([(1, 5), (1, 14), (4, 10), (7, 14), (7, 5)], [5, 23], [2]),  # M (0 - 2)
        ([(5, 2), (10, 2), (10, 8), (5, 8)], None, []),  # 0 - along edge - 0
        ([(5, 2), (10, 2), (10, 8), (15, 8)], None, []),  # 0 - along edge - 1
        ([(5, 2), (10, 2), (10, 14), (14, 14)], None, []),  # 0 - along edge - 2
        ([(5, 2), (10, 2), (10, 10), (8, 10), (8, 8)], None, []),
        ([(10, 2), (10, 5), (5, 5)], None, []),  # start at edge
        ([(10, 2), (10, 5), (5, 5), (5, 15)], None, []),  # start at edge - 0 - 2
        ([(10, 2), (10, 5), (5, 5), (5, 15), (15, 15)], None, [2]),
        ([(5, 5), (5, 10), (7, 10), (7, 15), (10, 15), (10, 17), (15, 17)], None, [2]),
        ([(5.6, 7), (10, 10.3), (15, 10.3)], [5.25], []),  # 0 to 3, 2 is below thresh
        ([(5.6, 7), (10, 10.3), (14.4, 7)], [5.5], []),  # 0 to 1, 2&3 are below thresh
        ([(0, 5), (18, 5)], [10], []),  # begins at model edge
        ([(0, 5), (0, 15), (5, 15)], [5], []),  # begins tangent to the model edge
        ([(0, 5), (18, 5), (18, 12)], [10, 23], [1]),  # begins at model edge, 2
        ([(0, 5), (0, 15), (15, 15)], [5, 20], [2]),  # begins tangent to model edge, 2
    ],
)
def test_embed_linear_object(grid2d, geometry, lines_s1d, embedded_in, reverse):
    if reverse:
        geometry = geometry[::-1]
        embedded_in = embedded_in[::-1]
        if lines_s1d is not None:
            lines_s1d = shapely.length(shapely.linestrings(geometry)) - lines_s1d[::-1]

    linear_objects = LinearObjects(
        id=[0],
        calculation_type=EMBEDDED,
        the_geom=[shapely.linestrings(geometry)],
    )

    nodes, actual_lines_s1d, line_ds1d_half = embed_linear_objects(
        linear_objects,
        grid2d.cell_tree,
        embedded_cutoff_threshold=1,
        embedded_node_id_counter=count(2),
    )

    assert_array_equal(nodes.embedded_in, embedded_in)
    if lines_s1d is not None:
        assert_almost_equal(actual_lines_s1d, lines_s1d)


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize(
    "geometry",
    [
        [(-1, 5), (18, 5)],  # begins outside of model
        [(5, 5), (-5, 5), (-5, 15), (5, 15)],  # begins in the model, but goes outside
    ],
)
def test_embed_linear_object_outside_raise(grid2d, geometry, reverse):
    if reverse:
        geometry = geometry[::-1]

    linear_objects = LinearObjects(
        id=[0],
        calculation_type=EMBEDDED,
        the_geom=[shapely.linestrings(geometry)],
    )

    with pytest.raises(
        SchematisationError,
        match=r"LinearObjects \[0\] are not completely inside the 2D cell.",
    ):
        embed_linear_objects(
            linear_objects,
            grid2d.cell_tree,
            embedded_cutoff_threshold=1,
            embedded_node_id_counter=count(2),
        )
