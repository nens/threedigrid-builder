from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import NodeType

import pytest


@pytest.fixture
def nodes():
    return Nodes(id=[1, 2, 3], coordinates=[(1, 10), (2, 12), (3, 13)])


def test_get_extent_1d(nodes):
    nodes = Nodes(
        id=range(3),
        coordinates=[(1, 10), (2, 12), (3, 13)],
        node_type=[
            NodeType.NODE_1D_STORAGE,
            NodeType.NODE_1D_NO_STORAGE,
            NodeType.NODE_2D_OPEN_WATER,
        ],
    )
    assert nodes.get_extent_1d() == (1, 10, 2, 12)


def test_get_extent_2d(nodes):
    nodes = Nodes(
        id=range(3),
        bounds=[
            (1, 10, 2, 11),
            (2, 12, 3, 13),
            (3, 13, 4, 14),
        ],
        node_type=[
            NodeType.NODE_2D_OPEN_WATER,
            NodeType.NODE_2D_OPEN_WATER,
            NodeType.NODE_1D_STORAGE,
        ],
    )
    assert nodes.get_extent_2d() == (1, 10, 3, 13)


def test_get_extent_none(nodes):
    nodes = Nodes(id=range(3))
    assert nodes.get_extent_1d() is None
    assert nodes.get_extent_2d() is None
