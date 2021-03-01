from numpy.testing import assert_array_equal
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import Grid

import numpy as np
import pygeos
import pytest
import itertools


@pytest.fixture
def connection_nodes():
    return ConnectionNodes(
        the_geom=pygeos.points([(0, 0), (10, 0)]),
        id=np.array([1, 3]),
        code=np.array(["one", "two"]),
    )


def test_get_grid(connection_nodes):
    counter = itertools.count(start=2)

    grid = Grid.from_connection_nodes(connection_nodes, counter)

    assert_array_equal(grid.nodes.id, [2, 3])
    assert next(counter) == 4
    assert_array_equal(grid.nodes.coordinates, [(0, 0), (10, 0)])
    assert_array_equal(grid.nodes.content_pk, [1, 3])
