import numpy as np
import pygeos
import pytest
from numpy.testing import assert_almost_equal, assert_equal

from threedigrid_builder.base import Endpoints, Lines, Nodes
from threedigrid_builder.constants import ContentType

CH = ContentType.TYPE_V2_CHANNEL
CN = ContentType.TYPE_V2_CONNECTION_NODES


@pytest.fixture
def nodes():
    return Nodes(
        id=[1, 2, 3], coordinates=[(1, 1), (2, 2), (3, 3)], content_type=[CN, CN, -9999]
    )


@pytest.fixture
def lines():
    return Lines(
        id=[0, 1, 2],
        line=[(1, 2), (2, 3), (1, 3)],
        line_geometries=[None, pygeos.linestrings([(5, 5), (6, 6)]), None],
        ds1d=[np.nan, 16.0, 2.0],
        content_type=[CN, CH, CH],
    )


def test_for_connection_nodes(nodes: Nodes, lines: Lines):
    endpoints = Endpoints.for_connection_nodes(nodes, lines)
    assert isinstance(endpoints, Endpoints)
    assert_equal(endpoints.node_id, [1, 1, 2, 2])
    assert_equal(endpoints.line_idx, [0, 2, 0, 1])
    assert_equal(endpoints.is_start, [1, 1, 0, 1])


def test_for_connection_nodes_with_line_type_filter(nodes: Nodes, lines: Lines):
    endpoints = Endpoints.for_connection_nodes(nodes, lines, line_types=[CH])
    assert isinstance(endpoints, Endpoints)
    assert_equal(endpoints.node_id, [1, 2])
    assert_equal(endpoints.line_idx, [2, 1])
    assert_equal(endpoints.is_start, [1, 1])


def test_getattr(lines: Lines):
    endpoints = Endpoints(lines=lines, id=[0, 1, 2], line_idx=[1, 2, 1])

    assert_equal(endpoints.ds1d, lines.ds1d[[1, 2, 1]])


@pytest.mark.parametrize(
    "ufunc,expected",
    [
        (np.fmax, [16.0, 16.0]),
        (np.fmin, [2.0, 16.0]),
        (np.add, [18.0, np.nan]),
    ],
)
def test_reduce_per_node(lines: Lines, ufunc, expected):
    endpoints = Endpoints(
        lines=lines, id=range(4), node_id=[1, 1, 5, 5], line_idx=[1, 2, 1, 0]
    )
    node_ids, maximums = endpoints.reduce_per_node(ufunc, endpoints.ds1d)

    assert_equal(node_ids, [1, 5])
    assert_almost_equal(maximums, expected)


def test_invert_level(lines: Lines):
    lines.invert_level_start_point[:] = [1, 2, 3]
    lines.invert_level_end_point[:] = [4, 5, 6]

    endpoints = Endpoints(
        lines=lines, id=range(3), line_idx=[0, 2, 1], is_start=[False, True, False]
    )

    assert_almost_equal(endpoints.invert_level, [4, 3, 5])