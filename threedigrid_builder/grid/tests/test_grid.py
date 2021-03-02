from numpy.testing import assert_array_equal
from numpy.testing import assert_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import CrossSectionLocations
from threedigrid_builder.grid import Grid

import itertools
import numpy as np
import pygeos
import pytest


@pytest.fixture
def connection_nodes():
    return ConnectionNodes(
        the_geom=pygeos.points([(0, 0), (10, 0)]),
        id=np.array([1, 3]),
        code=np.array(["one", "two"]),
    )


@pytest.fixture
def channel_grid():
    nodes_cn = Nodes(
        id=[0, 1, 2, 3],
        content_type=ContentType.TYPE_V2_CONNECTION_NODES,
        content_pk=[22, 24, 44, 66],
    )
    nodes_ch = Nodes(
        id=[4, 5],
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=[52, 52],
    )
    lines = Lines(
        id=[0, 1, 2, 3],
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=[51, 52, 52, 52],
        line=[(0, 1), (2, 4), (4, 5), (5, 3)],
        line_geometries=[
            pygeos.linestrings([(0, 0), (0, 10)]),
            pygeos.linestrings([(0, 0), (10, 0)]),
            pygeos.linestrings([(10, 0), (20, 0)]),
            pygeos.linestrings([(20, 0), (30, 0)]),
        ],
    )
    return Grid(nodes_cn + nodes_ch, lines)


@pytest.fixture
def cross_section_locations():
    # this fixture is closely related to channel_grid
    return CrossSectionLocations(
        id=[124, 125, 126],
        the_geom=pygeos.points([(0, 5), (12, 0), (30, 0)]),
        definition_id=[6, 7, 8],
        channel_id=[51, 52, 52],
    )


def test_from_connection_nodes(connection_nodes):
    counter = itertools.count(start=2)

    grid = Grid.from_connection_nodes(connection_nodes, counter)

    assert_array_equal(grid.nodes.id, [2, 3])
    assert next(counter) == 4
    assert_array_equal(grid.nodes.coordinates, [(0, 0), (10, 0)])
    assert_array_equal(grid.nodes.content_pk, [1, 3])


def test_set_cross_section_weights(channel_grid, cross_section_locations):
    # The test set has the following situation:

    # 1 channel with no added nodes, the cs location at the velocity point:
    # O --U-- O
    # 1 channel with 2 added nodes, with 2 cs locations:
    # O --X-- O U-X-- O --X-- U
    channel_grid.set_cross_section_weights(cross_section_locations)

    assert_equal(channel_grid.lines.cross1, [6, 7, 7, 7])
    assert_equal(channel_grid.lines.cross2, [6, 7, 8, 8])
    assert_equal(channel_grid.lines.cross_weights, [1.0, 1.0, 10 / 13, 5 / 18])
