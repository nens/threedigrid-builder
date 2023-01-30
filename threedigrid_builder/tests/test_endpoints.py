import numpy as np
import pytest
import shapely
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
        line_geometries=[None, shapely.linestrings([(5, 5), (6, 6)]), None],
        ds1d=[np.nan, 16.0, 2.0],
        ds1d_half=[np.nan, 8.0, 1.1],
        content_type=[CN, CH, CH],
    )


def test_for_connection_nodes(nodes: Nodes, lines: Lines):
    endpoints = Endpoints.for_connection_nodes(nodes, lines)
    endpoints.reorder_by("node_id")
    assert isinstance(endpoints, Endpoints)
    assert_equal(endpoints.node_id, [1, 1, 2, 2])
    assert_equal(endpoints.line_idx, [0, 2, 0, 1])
    assert_equal(endpoints.is_start, [1, 1, 0, 1])


def test_for_connection_nodes_with_line_mask(nodes: Nodes, lines: Lines):
    endpoints = Endpoints.for_connection_nodes(
        nodes, lines, line_mask=lines.content_type == CH
    )
    endpoints.reorder_by("node_id")
    assert isinstance(endpoints, Endpoints)
    assert_equal(endpoints.node_id, [1, 2])
    assert_equal(endpoints.line_idx, [2, 1])
    assert_equal(endpoints.is_start, [1, 1])


def test_for_connection_nodes_with_node_mask(nodes: Nodes, lines: Lines):
    endpoints = Endpoints.for_connection_nodes(nodes, lines, node_mask=nodes.id != 2)
    assert isinstance(endpoints, Endpoints)
    assert_equal(endpoints.node_id, [1, 1])
    assert_equal(endpoints.line_idx, [0, 2])
    assert_equal(endpoints.is_start, [1, 1])


def test_getattr(lines: Lines):
    endpoints = Endpoints(lines=lines, id=[0, 1, 2], line_idx=[1, 2, 1])

    assert_equal(endpoints.ds1d, lines.ds1d[[1, 2, 1]])


@pytest.mark.parametrize(
    "ufunc,expected",
    [
        (np.fmax.reduceat, [16.0, 16.0]),
        (np.fmin.reduceat, [2.0, 16.0]),
        (np.add.reduceat, [18.0, np.nan]),
        (lambda x, y: y, [0.0, 2.0]),
    ],
)
def test_reduce_per_node(lines: Lines, ufunc, expected):
    endpoints = Endpoints(
        lines=lines, id=range(4), node_id=[1, 1, 5, 5], line_idx=[1, 2, 1, 0]
    )
    actual = endpoints.reduce_per_node(ufunc, endpoints.ds1d)

    assert_almost_equal(actual, expected)


def test_reduce_per_node_2dim(lines: Lines):
    endpoints = Endpoints(
        lines=lines, id=range(4), node_id=[1, 1, 5, 5], line_idx=[1, 2, 1, 0]
    )
    actual = endpoints.reduce_per_node(
        np.fmax.reduceat, np.array([[1.0, 4.0], [3.0, 3.0], [5.0, 6.0], [7.0, 8.0]])
    )

    assert_almost_equal(actual, np.array([[3.0, 4.0], [7.0, 8.0]]))


def test_reduce_per_node_empty_2dim(lines: Lines):
    endpoints = Endpoints(lines=lines, id=[])
    actual = endpoints.reduce_per_node(None, np.empty((0, 2)))

    assert actual.shape == (0, 2)


def test_reduce_per_node_empty(lines: Lines):
    endpoints = Endpoints(lines=lines, id=[])
    actual = endpoints.reduce_per_node(None, [])

    assert actual.shape == (0,)


@pytest.fixture
def endpoints_for_agg(lines):
    return Endpoints(
        lines=lines, id=range(4), node_id=[1, 1, 5, 5], line_idx=[1, 2, 0, 1]
    )


def test_reduce_per_node_id(endpoints_for_agg: Endpoints):
    actual = endpoints_for_agg.get_reduce_per_node_id()

    assert_equal(actual, [1, 5])


def test_reduce_per_node_id_empty(lines: Lines):
    endpoints = Endpoints(lines=lines, id=[])
    actual = endpoints.get_reduce_per_node_id()

    assert actual.shape == (0,)


def test_nanmin_per_node(endpoints_for_agg: Endpoints):
    actual = endpoints_for_agg.nanmin_per_node(endpoints_for_agg.ds1d)

    assert_almost_equal(actual, [2.0, 16.0])


def test_first_per_node(endpoints_for_agg: Endpoints):
    actual = endpoints_for_agg.first_per_node(endpoints_for_agg.ds1d)

    assert_almost_equal(actual, [16.0, np.nan])


def test_sum_per_node(endpoints_for_agg: Endpoints):
    actual = endpoints_for_agg.sum_per_node(endpoints_for_agg.ds1d)

    assert_almost_equal(actual, [18.0, np.nan])


def test_invert_level(lines: Lines):
    lines.invert_level_start_point[:] = [1, 2, 3]
    lines.invert_level_end_point[:] = [4, 5, 6]

    endpoints = Endpoints(
        lines=lines, id=range(3), line_idx=[0, 2, 1], is_start=[False, True, False]
    )

    assert_almost_equal(endpoints.invert_level, [4, 3, 5])


def test_ds1d_endpoint(lines: Lines):
    endpoints = Endpoints(
        lines=lines, id=[0, 1, 2], line_idx=[1, 2, 2], is_start=[True, True, False]
    )

    assert_almost_equal(endpoints.ds1d_endpoint, [8.0, 1.1, 0.9])
