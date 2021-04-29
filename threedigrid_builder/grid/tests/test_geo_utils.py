from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal
from threedigrid_builder.grid.geo_utils import segmentize, segmentize_and_substring2, segmentize_and_substring

import numpy as np
import pygeos
import pytest


@pytest.fixture
def random_lines():
    n_lines = 10  # approximate value
    n_coords = 100
    indices = np.sort(np.random.randint(0, n_lines, n_coords))
    # discard lines with fewer than 2 coordinates
    indices = indices[~np.isin(indices, np.where(np.bincount(indices) < 2)[0])]
    lines = pygeos.linestrings(np.random.random((len(indices), 2)), indices=indices)
    # discard nones
    lines = lines[~pygeos.is_missing(lines)]
    return lines


@pytest.mark.parametrize(
    "segment_size,expected_size,expected_nodes",
    [
        (3.0, 12 / 4, [(3, 0), (6, 0), (6, 3)]),
        (4.0, 12 / 3, [(4, 0), (6, 2)]),
        (5.0, 12 / 2, [(6, 0)]),
        (6.0, 12 / 2, [(6, 0)]),
        (8.0, 12 / 2, [(6, 0)]),
        (9.0, 12, np.empty((0, 2), dtype=float)),
        (100.0, 12, np.empty((0, 2), dtype=float)),
    ],
)
def test_segmentize_one_line(segment_size, expected_size, expected_nodes):
    line = pygeos.linestrings([[(0, 0), (6, 0), (6, 6)]])
    nodes, node_idx, actual_size, segments, segment_idx = segmentize_and_substring(line, segment_size)

    assert_almost_equal(expected_nodes, pygeos.get_coordinates(nodes), decimal=7)
    # assert_array_equal(expected_lines, pygeos.get_coordinates(lines))
    assert_array_equal(segment_idx, 0)
    assert_array_equal(node_idx, 0)
    assert_almost_equal(actual_size, expected_size, decimal=7)


def test_segmentize_many(random_lines):
    nodes, node_idx, actual_size, segments, segment_idx = segmentize_and_substring(random_lines, 1.0)

    # the node coordinates should match line start & end coordinates
    idx_prev = -1
    line_cur = -1
    for _idx, _node in zip(node_idx, nodes):
        if _idx != idx_prev:
            line_cur += 1
            idx_prev = _idx
        node_coord = [pygeos.get_x(_node), pygeos.get_y(_node)]
        line_before = pygeos.get_coordinates(segments[line_cur])
        line_after = pygeos.get_coordinates(segments[line_cur + 1])
        assert_almost_equal(node_coord, line_before[-1])
        assert_almost_equal(node_coord, line_after[0])
        line_cur += 1

    # the lengths of each line segments should equal the reported size
    for _idx, _line in zip(segment_idx, segments):
        assert pygeos.length(_line) == pytest.approx(actual_size[_idx])

    # the number of segments times the reported size should equal the input line's size
    actual_size_total = (
        np.bincount(segment_idx, minlength=len(random_lines)) * actual_size
    )
    assert_almost_equal(actual_size_total, pygeos.length(random_lines), decimal=7)
