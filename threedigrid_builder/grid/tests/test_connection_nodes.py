from numpy.testing import assert_array_equal
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import NodeType
from threedigrid_builder.grid import ConnectionNodes

import itertools
import numpy as np
import pygeos
import pytest


@pytest.fixture
def connection_nodes():
    return ConnectionNodes(
        the_geom=pygeos.points([(0, 0), (10, 0), (10, 20), (10, 30)]),
        id=np.array([1, 3, 4, 9]),
        code=np.array(["one", "two", "", ""]),
        storage_area=np.array([-2, 0, 15, np.nan]),
        calculation_type=np.array([1, 2, 5, 2]),
        bottom_level=np.array([2.1, 5.2, np.nan, np.nan]),
    )


def test_get_nodes(connection_nodes):
    counter = itertools.count(start=2)

    nodes = connection_nodes.get_nodes(counter)
    assert isinstance(nodes, Nodes)

    assert_array_equal(nodes.id, [2, 3, 4, 5])
    assert next(counter) == 6
    assert_array_equal(nodes.coordinates, [(0, 0), (10, 0), (10, 20), (10, 30)])
    assert_array_equal(nodes.content_pk, [1, 3, 4, 9])
    assert_array_equal(
        nodes.node_type,
        [
            NodeType.NODE_1D_NO_STORAGE,
            NodeType.NODE_1D_NO_STORAGE,
            NodeType.NODE_1D_STORAGE,
            NodeType.NODE_1D_NO_STORAGE,
        ],
    )
    assert_array_equal(nodes.calculation_type, connection_nodes.calculation_type)
    assert_array_equal(nodes.dmax, connection_nodes.bottom_level)
