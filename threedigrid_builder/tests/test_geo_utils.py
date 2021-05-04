from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal
from threedigrid_builder.grid import geo_utils

import numpy as np
import pygeos
import pytest


@pytest.fixture(scope="session")
def random_lines():
    n_lines = 10  # approximate value
    n_coords = 100
    np.random.seed(0)
    indices = np.sort(np.random.randint(0, n_lines, n_coords))
    # discard lines with fewer than 2 coordinates
    indices = indices[~np.isin(indices, np.where(np.bincount(indices) < 2)[0])]
    coords = np.around(np.random.random((len(indices), 2)), decimals=6)
    lines = pygeos.linestrings(coords, indices=indices)
    # discard nones
    lines = lines[~pygeos.is_missing(lines)]
    return lines


@pytest.fixture
def two_lines():
    return np.array(
        [
            pygeos.linestrings([(0, 10), (10, 10)]),
            pygeos.linestrings([(0, 0), (6, 0), (6, 6)]),
        ]
    )


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
def test_segmentize_one(segment_size, expected_size, expected_nodes):
    line = pygeos.linestrings([[(0, 0), (6, 0), (6, 6)]])
    nodes, node_idx, s = geo_utils.segmentize(line, segment_size)

    assert_almost_equal(expected_nodes, pygeos.get_coordinates(nodes), decimal=7)
    assert_array_equal(node_idx, 0)
    assert_almost_equal(s, expected_size * (np.arange(len(expected_nodes)) + 1))


def test_segmentize_two(two_lines):
    nodes, node_idx, s = geo_utils.segmentize(two_lines, 5.0)

    assert_almost_equal(
        [[5.0, 10.0], [6.0, 0.0]], pygeos.get_coordinates(nodes), decimal=7
    )
    assert_array_equal(node_idx, [0, 1])
    assert_almost_equal(s, [5.0, 6.0])


def test_segmentize_many(random_lines):
    d = 1.0
    nodes, node_idx, s = geo_utils.segmentize(random_lines, d)

    assert nodes.shape == node_idx.shape

    # check the number of nodes per line
    n_segments = np.round(pygeos.length(random_lines) / d).astype(int)
    n_segments[n_segments < 1] = 1
    assert_array_equal(np.repeat(np.arange(len(n_segments)), n_segments - 1), node_idx)


@pytest.mark.parametrize(
    "start,end,expected_coords",
    [
        (0.0, 12.0, [(0.0, 0.0), (6.0, 0.0), (6.0, 6.0)]),
        (1.45, 12.0, [(1.45, 0.0), (6.0, 0.0), (6.0, 6.0)]),
        (6.0, 12.0, [(6.0, 0.0), (6.0, 6.0)]),
        (6.7, 12.0, [(6.0, 0.7), (6.0, 6.0)]),
        (0.0, 11.5, [(0.0, 0.0), (6.0, 0.0), (6.0, 5.5)]),
        (0.0, 6.0, [(0.0, 0.0), (6.0, 0.0)]),
        (0.0, 3.5, [(0.0, 0.0), (3.5, 0.0)]),
        (1.0, 9.0, [(1.0, 0.0), (6.0, 0.0), (6.0, 3.0)]),
        (6.0 + 1e-7, 12.0, [(6.0 + 1e-7, 0.0), (6.0, 6.0)]),
        (6.0 + 1e-8, 12.0, [(6.0, 0.0), (6.0, 6.0)]),
        (6.0, 12.0 + 1e-15, [(6.0, 0.0), (6.0, 6.0)]),
    ],
)
def test_line_substring_one(start, end, expected_coords):
    line = pygeos.linestrings([[(0, 0), (6, 0), (6, 6)]])
    segments = geo_utils.line_substring(line, start, end)

    assert_almost_equal(expected_coords, pygeos.get_coordinates(segments), decimal=7)


def test_line_substring_two(two_lines):
    segments = geo_utils.line_substring(two_lines, [3.0, 4.0], [9.0, 10.0])

    coords, index = pygeos.get_coordinates(segments, return_index=True)

    assert_almost_equal(
        coords, [(3.0, 10.0), (9.0, 10.0), (4.0, 0), (6.0, 0.0), (6.0, 4.0)]
    )
    assert_almost_equal(index, [0, 0, 1, 1, 1])


def test_line_substring_many(random_lines):
    nodes, node_idx, _ = geo_utils.segmentize(random_lines, 1.0)
    segment_counts = np.bincount(node_idx, minlength=len(random_lines)) + 1
    start, end, segment_idx = geo_utils.segment_start_end(random_lines, segment_counts)
    segments = geo_utils.line_substring(random_lines, start, end, segment_idx)

    # the length of each line segment should equal the line length / number of segments
    expected_sizes = pygeos.length(random_lines) / segment_counts
    for segment, idx in zip(segments, segment_idx):
        assert pygeos.length(segment) == pytest.approx(expected_sizes[idx])

    # check the start (resp. end) coordinate of the first (resp. last) segment per line
    segment_start_idx, segment_last_idx = geo_utils.counts_to_ranges(segment_counts)
    segment_last_idx -= 1
    for start, last, line in zip(segment_start_idx, segment_last_idx, random_lines):
        line_coords = pygeos.get_coordinates(line)
        assert_almost_equal(pygeos.get_coordinates(segments[start])[0], line_coords[0])
        assert_almost_equal(pygeos.get_coordinates(segments[last])[-1], line_coords[-1])
    # check the remaining segment start/end with the nodes
    node_idx = 0
    for idx, segment in enumerate(segments):
        segment_coords = pygeos.get_coordinates(segment)
        if idx not in segment_start_idx:
            assert_almost_equal(
                segment_coords[0], pygeos.get_coordinates(nodes[node_idx])[0]
            )
            node_idx += 1
        if idx not in segment_last_idx:
            assert_almost_equal(
                segment_coords[-1], pygeos.get_coordinates(nodes[node_idx])[0]
            )
    # check the segment mid-coordinates with the linestrings
    node_idx = 0
    for line_idx, segment in zip(segment_idx, segments):
        segment_coords = pygeos.get_coordinates(segment)
        line_coords = pygeos.get_coordinates(random_lines[line_idx])
        for coord in segment_coords[1:-1]:
            # segment mid coord should match 1 line coord
            assert np.isclose(coord, line_coords).all(axis=1).any()
