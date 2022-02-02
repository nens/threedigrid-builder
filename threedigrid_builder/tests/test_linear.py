from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal
from pygeos.testing import assert_geometries_equal
from threedigrid_builder.base import array_of
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import linear

import itertools
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


@pytest.fixture
def connection_nodes():
    # Used to map connection_node_start/end_id to an index (sequence id)
    return ConnectionNodes(
        id=np.array([21, 25, 33, 42]),
        the_geom=pygeos.points([(0, 21), (1, 25), (2, 33), (3, 42)]),
        drain_level=[0.0, 20.0, 30.0, 10.0],
    )


class LinearObject:
    id: int
    the_geom: pygeos.Geometry
    calculation_type: CalculationType
    connection_node_start_id: int
    connection_node_end_id: int
    dist_calc_points: float
    cross_section_definition_id: int
    invert_level_start_point: float
    invert_level_end_point: float
    discharge_coefficient_positive: float
    discharge_coefficient_negative: float
    display_name: str
    zoom_category: int
    material: int
    sewerage_type: int


@array_of(LinearObject)
class LinearObjects(linear.BaseLinear):
    content_type = ContentType.TYPE_V2_CONNECTION_NODES  # just pick one for the test


@pytest.fixture
def one_linear_object():
    return LinearObjects(
        the_geom=pygeos.linestrings([[(0, 0), (6, 0), (6, 6)]]),
        dist_calc_points=np.array([5.0]),
        id=np.array([1]),
        connection_node_start_id=np.array([21]),
        connection_node_end_id=np.array([42]),
        calculation_type=np.array([2]),
        cross_section_definition_id=[17],
    )


@pytest.fixture
def two_linear_objects():
    return LinearObjects(
        the_geom=pygeos.linestrings(
            [[(0, 0), (10, 0), (10, 10)], [(0, 0), (0, 100), (100, 100)]]
        ),
        dist_calc_points=[5.0, np.nan],
        id=[1, 2],
        connection_node_start_id=[21, 25],
        connection_node_end_id=[42, 33],
        calculation_type=[2, 1],
        cross_section_definition_id=[17, 18],
        invert_level_start_point=[1.0, 2.0],
        invert_level_end_point=[3.0, 4.0],
        discharge_coefficient_positive=[0.8, 0.5],
        discharge_coefficient_negative=[0.5, 0.8],
        material=[-9999, 1],
        sewerage_type=[2, 1],
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
    nodes, node_idx, s = linear.segmentize(line, segment_size)

    assert_almost_equal(expected_nodes, pygeos.get_coordinates(nodes), decimal=7)
    assert_array_equal(node_idx, 0)
    assert_almost_equal(s, expected_size * (np.arange(len(expected_nodes)) + 1))


def test_segmentize_two(two_lines):
    nodes, node_idx, s = linear.segmentize(two_lines, 5.0)

    assert_almost_equal(
        [[5.0, 10.0], [6.0, 0.0]], pygeos.get_coordinates(nodes), decimal=7
    )
    assert_array_equal(node_idx, [0, 1])
    assert_almost_equal(s, [5.0, 6.0])


def test_segmentize_many(random_lines):
    d = 1.0
    nodes, node_idx, s = linear.segmentize(random_lines, d)

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
    segments = linear.line_substring(line, start, end)

    assert_almost_equal(expected_coords, pygeos.get_coordinates(segments), decimal=7)


def test_line_substring_two(two_lines):
    segments = linear.line_substring(two_lines, [3.0, 4.0], [9.0, 10.0])

    coords, index = pygeos.get_coordinates(segments, return_index=True)

    assert_almost_equal(
        coords, [(3.0, 10.0), (9.0, 10.0), (4.0, 0), (6.0, 0.0), (6.0, 4.0)]
    )
    assert_almost_equal(index, [0, 0, 1, 1, 1])


def test_line_substring_many(random_lines):
    nodes, node_idx, dist_to_start = linear.segmentize(random_lines, 1.0)
    segment_counts = np.bincount(node_idx, minlength=len(random_lines)) + 1
    start, end, segment_idx = linear.segment_start_end(
        random_lines, segment_counts, dist_to_start
    )
    segments = linear.line_substring(random_lines, start, end, segment_idx)

    # the length of each line segment should equal the line length / number of segments
    expected_sizes = pygeos.length(random_lines) / segment_counts
    for segment, idx in zip(segments, segment_idx):
        assert pygeos.length(segment) == pytest.approx(expected_sizes[idx])

    # check the start (resp. end) coordinate of the first (resp. last) segment per line
    segment_start_idx, segment_last_idx = linear.counts_to_ranges(segment_counts)
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


@pytest.mark.parametrize(
    "dist,expected",
    [
        (3.0, [(3, 0), (6, 0), (6, 3)]),  # 12 / 3  = 4 segments
        (4.0, [(4, 0), (6, 2)]),  # 12 / 4  = 3 segments
        (5.0, [(6, 0)]),  # 12 / 5  = 2.4 -> 2 segments
        (6.0, [(6, 0)]),  # 12 / 6  = 2 segments
        (8.0, [(6, 0)]),  # 12 / 8  = 1.5 -> 2
        (9.0, np.empty((0, 2), dtype=float)),  # 12 / 9  = 1.33 -> 1
        (100.0, np.empty((0, 2), dtype=float)),
        (0.0, np.empty((0, 2), dtype=float)),
        (-9999.0, np.empty((0, 2), dtype=float)),
        (np.nan, np.empty((0, 2), dtype=float)),
    ],
)
def test_interpolate_nodes_one_linear_object(dist, expected, one_linear_object):
    one_linear_object.dist_calc_points[0] = dist
    nodes = one_linear_object.interpolate_nodes(
        itertools.count(start=2), global_dist_calc_points=74.0
    )

    assert_array_equal(nodes.id, range(2, 2 + len(expected)))
    assert_almost_equal(nodes.coordinates, expected, decimal=7)
    assert_array_equal(nodes.content_type, ContentType.TYPE_V2_CONNECTION_NODES)
    assert_array_equal(nodes.content_pk, 1)
    assert_array_equal(nodes.node_type, NodeType.NODE_1D_NO_STORAGE)
    assert_array_equal(nodes.calculation_type, 2)

    expected_size = 12.0 / (len(expected) + 1)
    assert_array_equal(nodes.s1d, np.arange(1, len(expected) + 1) * expected_size)


@pytest.mark.parametrize("global_value", [0.0, -9999.0, np.nan, None])
def test_interpolate_nodes_global_none(one_linear_object, global_value):
    expected = np.empty((0, 2), dtype=float)

    one_linear_object.dist_calc_points[0] = np.nan
    nodes = one_linear_object.interpolate_nodes(
        itertools.count(start=2), global_dist_calc_points=global_value
    )

    assert_array_equal(nodes.id, range(2, 2 + len(expected)))
    assert_almost_equal(nodes.coordinates, expected, decimal=7)
    assert_array_equal(nodes.content_type, ContentType.TYPE_V2_CONNECTION_NODES)
    assert_array_equal(nodes.content_pk, 1)
    assert_array_equal(nodes.node_type, NodeType.NODE_1D_NO_STORAGE)
    assert_array_equal(nodes.calculation_type, 2)

    expected_size = 12.0 / (len(expected) + 1)
    assert_array_equal(nodes.s1d, np.arange(1, len(expected) + 1) * expected_size)


def test_interpolate_nodes_two_linear_objects(two_linear_objects):
    nodes = two_linear_objects.interpolate_nodes(
        itertools.count(start=2), global_dist_calc_points=50.0
    )

    expected_points = [(5, 0), (10, 0), (10, 5), (0, 50), (0, 100), (50, 100)]

    assert_array_equal(nodes.id, range(2, 8))
    assert_array_equal(nodes.coordinates, expected_points)
    assert_array_equal(nodes.content_type, ContentType.TYPE_V2_CONNECTION_NODES)
    assert_array_equal(nodes.content_pk, [1, 1, 1, 2, 2, 2])
    assert_array_equal(nodes.node_type, NodeType.NODE_1D_NO_STORAGE)
    assert_array_equal(nodes.calculation_type, [2, 2, 2, 1, 1, 1])


def test_interpolate_nodes_skips_embedded(two_linear_objects):
    two_linear_objects.calculation_type[0] = CalculationType.EMBEDDED
    nodes = two_linear_objects.interpolate_nodes(
        itertools.count(start=2), global_dist_calc_points=50.0
    )

    assert_array_equal(nodes.coordinates, [(0, 50), (0, 100), (50, 100)])
    assert_array_equal(nodes.content_pk, 2)
    assert_array_equal(nodes.calculation_type, 1)


def test_get_lines(connection_nodes, two_linear_objects):
    nodes = Nodes(
        id=[10, 11, 12],
        content_pk=[1, 2, 2],
        s1d=[10.0, 50.0, 150.0],  # segments are not precisely the same size!
    )

    lines = two_linear_objects.get_lines(
        connection_nodes,
        nodes,
        itertools.count(start=0),
        connection_node_offset=100,
    )

    expected_line = [(100, 10), (10, 103), (101, 11), (11, 12), (12, 102)]
    expected_centers = [5.0, 15.0, 25.0, 100.0, 175.0]
    expected_sizes = [10.0, 10.0, 50.0, 100.0, 50.0]

    assert_array_equal(lines.id, range(5))
    assert_array_equal(lines.line, expected_line)
    assert_array_equal(lines.content_type, ContentType.TYPE_V2_CONNECTION_NODES)
    assert_array_equal(lines.content_pk, [1, 1, 2, 2, 2])
    assert_array_equal(lines.kcu, [2, 2, 1, 1, 1])
    assert_array_equal(lines.cross_id1, [17, 17, 18, 18, 18])
    assert_array_equal(lines.cross_id2, [17, 17, 18, 18, 18])
    assert_array_equal(lines.cross_weight, 1.0)
    assert_array_equal(lines.discharge_coefficient_positive, [0.8, 1.0, 0.5, 1.0, 1.0])
    assert_array_equal(lines.discharge_coefficient_negative, [1.0, 0.5, 1.0, 1.0, 0.8])
    assert_almost_equal(lines.s1d, expected_centers)
    assert_almost_equal(lines.ds1d, expected_sizes)
    assert_almost_equal(pygeos.length(lines.line_geometries), expected_sizes)
    assert_almost_equal(lines.invert_level_start_point, [1.0, 2.0, 2.0, 2.5, 3.5])
    assert_almost_equal(lines.invert_level_end_point, [2.0, 3.0, 2.5, 3.5, 4])
    assert_almost_equal(lines.dpumax, [2, 3, 2.5, 3.5, 4])
    assert_almost_equal(lines.material, [-9999, -9999, 1, 1, 1])
    assert_almost_equal(lines.sewerage_type, [2, 2, 1, 1, 1])


def test_get_lines_embedded_mode(connection_nodes, two_linear_objects):
    nodes = Nodes(
        id=[10, 11, 12],
        content_pk=[1, 2, 2],
        s1d=[0.1] * 3,  # anything that isn't NaN
        embedded_in=[5, 4, 1],
    )

    two_linear_objects.calculation_type[:] = CalculationType.EMBEDDED
    lines = two_linear_objects.get_lines(
        connection_nodes,
        nodes,
        itertools.count(start=0),
        connection_node_offset=100,
        embedded_mode=True,
    )

    assert_array_equal(lines.line, [(100, 5), (5, 103), (101, 4), (4, 1), (1, 102)])


def test_get_lines_non_embedded_mode_skips(connection_nodes, two_linear_objects):
    # embedded linear objects are skipped
    two_linear_objects.calculation_type[1] = CalculationType.EMBEDDED
    lines = two_linear_objects.get_lines(
        connection_nodes,
        Nodes(id=[]),
        itertools.count(start=0),
        connection_node_offset=100,
        embedded_mode=False,
    )
    assert_array_equal(lines.line, [(100, 103)])


def test_get_lines_embedded_mode_skips(connection_nodes, two_linear_objects):
    # non-embedded linear objects are skipped
    two_linear_objects.calculation_type[1] = CalculationType.EMBEDDED
    lines = two_linear_objects.get_lines(
        connection_nodes,
        Nodes(id=[]),
        itertools.count(start=0),
        connection_node_offset=100,
        embedded_mode=True,
    )
    assert_array_equal(lines.line, [(101, 102)])


@pytest.mark.parametrize(
    "linear_object_idx,expected",
    [
        ([], [(0, 3)]),
        ([1], [(0, 4), (4, 3)]),
        ([1, 1, 1], [(0, 4), (4, 5), (5, 6), (6, 3)]),
    ],
)
def test_get_lines_one_linear_object(
    linear_object_idx, expected, connection_nodes, one_linear_object
):
    nodes = Nodes(
        id=range(4, 4 + len(linear_object_idx)),
        content_pk=linear_object_idx,
        s1d=[0.1] * (len(linear_object_idx)),  # some number, doesn't matter
    )
    lines = one_linear_object.get_lines(
        connection_nodes, nodes, itertools.count(start=0)
    )

    assert_array_equal(lines.line, expected)


@pytest.mark.parametrize(
    "linear_object_idx,expected",
    [
        ([], [(0, 3), (1, 2)]),
        ([1], [(0, 4), (4, 3), (1, 2)]),
        ([1, 1, 1], [(0, 4), (4, 5), (5, 6), (6, 3), (1, 2)]),
        ([2], [(0, 3), (1, 4), (4, 2)]),
        ([2, 2, 2], [(0, 3), (1, 4), (4, 5), (5, 6), (6, 2)]),
        ([1, 2, 2], [(0, 4), (4, 3), (1, 5), (5, 6), (6, 2)]),
        ([1, 1, 2], [(0, 4), (4, 5), (5, 3), (1, 6), (6, 2)]),
    ],
)
def test_get_lines_two_linear_objects(
    linear_object_idx, expected, connection_nodes, two_linear_objects
):
    nodes = Nodes(
        id=range(4, 4 + len(linear_object_idx)),
        content_pk=linear_object_idx,
        s1d=[0.1] * (len(linear_object_idx)),  # some number, doesn't matter
    )

    lines = two_linear_objects.get_lines(
        connection_nodes, nodes, itertools.count(start=0)
    )

    assert_array_equal(lines.line, expected)


@pytest.mark.parametrize(
    "geom,expected",
    [
        (None, [(0, 21), (3, 42)]),
        ([(0, 21), (3, 42)], [(0, 21), (3, 42)]),
        ([(0, 21), (1, 2), (3, 42)], [(0, 21), (1, 2), (3, 42)]),
        ([(0, 21), (3, 42)][::-1], [(0, 21), (3, 42)]),
        ([(0, 21), (1, 2), (3, 42)][::-1], [(0, 21), (1, 2), (3, 42)]),
        ([(1, 21), (3, 41)], [(0, 21), (1, 21), (3, 41), (3, 42)]),
        ([(3, 41), (1, 21)], [(0, 21), (1, 21), (3, 41), (3, 42)]),
    ],
)
def test_set_geometries(two_linear_objects, connection_nodes, geom, expected):
    two_linear_objects.the_geom[0] = pygeos.linestrings(geom) if geom else None
    two_linear_objects.set_geometries(connection_nodes)

    expected_geometry = pygeos.linestrings(expected)
    assert_geometries_equal(two_linear_objects.the_geom[0], expected_geometry)


def test_interpolate_nodes_no_geometries(two_linear_objects):
    two_linear_objects.the_geom[:] = None
    with pytest.raises(ValueError, match=".*encountered without a geometry."):
        two_linear_objects.interpolate_nodes(None, None)


@pytest.mark.parametrize(
    "ids,ds,expected",
    [
        ([1], [10.0], [5.0]),  # halfway the first object
        ([2], [100.0], [25.0]),  # halfway the second object
        ([1, 2], [10.0, 100.0], [5.0, 25.0]),
        ([2, 1], [100.0, 10.0], [25.0, 5.0]),
        ([1, 1, 1], [0.0, 5.0, 20.0], [0.0, 2.5, 10.0]),  # at 0%, 25%, 100%
        ([1, 2, 1], [0.0, 150.0, 20.0], [0.0, 27.5, 10.0]),  # mixed
    ],
)
def test_compute_drain_level(ids, ds, expected, two_linear_objects, connection_nodes):
    actual = two_linear_objects.compute_drain_level(ids, ds, connection_nodes)

    assert_almost_equal(actual, expected)


def test_compute_drain_level_with_one_nan(two_linear_objects, connection_nodes):
    connection_nodes.drain_level[:] = [0.0, np.nan, 1.0, np.nan]
    actual = two_linear_objects.compute_drain_level(
        [1, 2], [3.0, 4.0], connection_nodes
    )

    assert_almost_equal(actual, [0.0, 1.0])


def test_compute_drain_level_with_two_nan(two_linear_objects, connection_nodes):
    connection_nodes.drain_level[:] = [np.nan, 1.0, 1.0, np.nan]
    actual = two_linear_objects.compute_drain_level(
        [1, 2], [3.0, 4.0], connection_nodes
    )

    assert_almost_equal(actual, [np.nan, 1.0])
