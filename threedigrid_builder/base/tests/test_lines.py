from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes

import pygeos
import pytest


@pytest.fixture
def nodes():
    return Nodes(id=[1, 2, 3], coordinates=[(1, 1), (2, 2), (3, 3)])


@pytest.fixture
def lines():
    return Lines(
        id=[0, 1, 2],
        line=[(1, 2), (2, 3), (1, 3)],
        line_geometries=[None, pygeos.linestrings([(5, 5), (6, 6)]), None],
    )


def test_fix_line_geometries(nodes, lines):
    lines.fix_line_geometries(nodes)

    assert pygeos.equals(lines.line_geometries[0], pygeos.linestrings([(1, 1), (2, 2)]))
    assert pygeos.equals(lines.line_geometries[1], pygeos.linestrings([(5, 5), (6, 6)]))
    assert pygeos.equals(lines.line_geometries[2], pygeos.linestrings([(1, 1), (3, 3)]))
