import itertools
from unittest import mock

import numpy as np
import pytest
import shapely
from numpy.testing import assert_almost_equal, assert_array_equal
from shapely.testing import assert_geometries_equal

from threedigrid_builder.base import Array, Lines, Nodes, PointsOnLine
from threedigrid_builder.constants import CalculationType, ContentType, NodeType
from threedigrid_builder.grid import ConnectionNodes, linear


@pytest.fixture(scope="session")
def random_lines():
    n_lines = 10  # approximate value
    n_coords = 100
    np.random.seed(0)
    indices = np.sort(np.random.randint(0, n_lines, n_coords))
    # discard lines with fewer than 2 coordinates
    indices = indices[~np.isin(indices, np.where(np.bincount(indices) < 2)[0])]
    coords = np.around(np.random.random((len(indices), 2)), decimals=6)
    lines = shapely.linestrings(coords, indices=indices)
    # discard nones
    lines = lines[~shapely.is_missing(lines)]
    return lines


@pytest.fixture
def two_lines():
    return np.array(
        [
            shapely.linestrings([(0, 10), (10, 10)]),
            shapely.linestrings([(0, 0), (6, 0), (6, 6)]),
        ]
    )


@pytest.fixture
def connection_nodes():
    # Used to map connection_node_start/end_id to an index (sequence id)
    return ConnectionNodes(
        id=np.array([21, 25, 33, 42]),
        the_geom=shapely.points([(0, 21), (1, 25), (2, 33), (3, 42)]),
        drain_level=[0.0, 20.0, 30.0, 10.0],
    )


class LinearObject:
    id: int
    the_geom: shapely.Geometry
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


class LinearObjects(Array[LinearObject], linear.BaseLinear):
    content_type = ContentType.TYPE_V2_1D_LATERAL  # just pick one for the test


@pytest.fixture
def one_linear_object():
    return LinearObjects(
        the_geom=shapely.linestrings([[(0, 0), (6, 0), (6, 6)]]),
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
        the_geom=shapely.linestrings(
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
    "global_dist_calc_points,expected_dists",
    [
        (74.0, [5.0, 74.0]),
        (0.0, [5.0, np.inf]),
        (-9999.0, [5.0, np.inf]),
        (np.nan, [5.0, np.inf]),
        (None, [5.0, np.inf]),
    ],
)
def test_interpolate_nodes(global_dist_calc_points, expected_dists, two_linear_objects):
    dummy_points = PointsOnLine(
        linestrings=two_linear_objects.linestrings,
        id=[0, 1, 2],
        s1d=[1.0, 2.0, 3.0],
        linestring_idx=[0, 0, 1],
    )
    with mock.patch.object(two_linear_objects.__class__, "linestrings") as linestrings:
        linestrings.segmentize().linestring_idx = [0, 1]
        linestrings.segmentize().interpolate_points.return_value = dummy_points
        nodes = two_linear_objects.interpolate_nodes(
            itertools.count(start=2), global_dist_calc_points
        )

    (fixed_nodes,), _ = linestrings.segmentize.call_args
    assert isinstance(fixed_nodes, PointsOnLine)
    assert len(fixed_nodes) == 0

    (dists,), _ = linestrings.segmentize().interpolate_points.call_args
    assert_almost_equal(dists, expected_dists)

    assert_array_equal(nodes.id, range(2, 2 + len(dummy_points)))
    assert_almost_equal(
        nodes.coordinates, shapely.get_coordinates(dummy_points.the_geom)
    )
    assert_array_equal(nodes.content_type, ContentType.TYPE_V2_1D_LATERAL)
    assert_array_equal(nodes.content_pk, [1, 1, 2])
    assert_array_equal(nodes.node_type, NodeType.NODE_1D_NO_STORAGE)
    assert_array_equal(nodes.calculation_type, [2, 2, 1])
    assert_array_equal(nodes.s1d, dummy_points.s1d)


def test_interpolate_nodes_skips_embedded(two_linear_objects):
    two_linear_objects.calculation_type[0] = CalculationType.EMBEDDED

    with mock.patch.object(LinearObjects, "linestrings") as linestrings:
        linestrings.segmentize().linestring_idx = [0, 1]
        linestrings.segmentize().interpolate_points.return_value = PointsOnLine.empty(
            two_linear_objects.linestrings
        )
        two_linear_objects.interpolate_nodes(
            itertools.count(start=2), global_dist_calc_points=74.0
        )

    (dists,), _ = linestrings.segmentize().interpolate_points.call_args
    assert_almost_equal(dists, [np.inf, 74.0])


def test_interpolate_nodes_with_fixture(two_linear_objects):
    fixed_nodes = PointsOnLine(
        linestrings=two_linear_objects.linestrings,
        id=[0, 1],
        s1d=[30.0, 50.0],
        linestring_idx=[1, 1],
        content_pk=[1, 2],
        secondary_content_pk=[-9999, 4],
    )
    dummy_points = PointsOnLine(
        linestrings=two_linear_objects.linestrings,
        id=[0, 1],
        s1d=[10.0, 125.0],
        linestring_idx=[0, 1],
    )
    with mock.patch.object(LinearObjects, "linestrings") as linestrings:
        linestrings.segmentize().linestring_idx = np.array([0, 1, 1, 1])
        linestrings.segmentize().interpolate_points.return_value = dummy_points
        nodes = two_linear_objects.interpolate_nodes(
            itertools.count(start=2),
            global_dist_calc_points=74.0,
            fixed_nodes=fixed_nodes,
        )

    (arg,), _ = linestrings.segmentize.call_args
    assert isinstance(arg, PointsOnLine)
    assert_array_equal(arg.s1d, fixed_nodes.s1d)

    (dists,), _ = linestrings.segmentize().interpolate_points.call_args
    assert_almost_equal(dists, [5.0, 74.0, 74.0, 74.0])

    assert_array_equal(nodes.id, [2, 3, 4, 5])
    assert_almost_equal(nodes.coordinates, [(10, 0), (0, 30), (0, 50), (25, 100)])
    assert_array_equal(nodes.content_type, ContentType.TYPE_V2_1D_LATERAL)
    assert_array_equal(nodes.content_pk, [1, 2, 2, 2])
    assert_array_equal(nodes.calculation_type, [2, 1, 1, 1])
    assert_array_equal(nodes.s1d, [10.0, 30.0, 50.0, 125.0])
    assert_array_equal(
        nodes.breach_ids, [(-9999, -9999), (1, -9999), (2, 4), (-9999, -9999)]
    )


def test_interpolate_nodes_with_fixture_on_start_end(one_linear_object):
    fixed_nodes = PointsOnLine(
        linestrings=one_linear_object.linestrings,
        id=[0, 1],
        s1d=[0.0, 12.0],
        linestring_idx=[0, 0],
        content_pk=[1, 2],
        secondary_content_pk=[-9999, 4],
    )
    dummy_points = PointsOnLine.empty(one_linear_object.linestrings)
    with mock.patch.object(LinearObjects, "linestrings") as linestrings:
        linestrings.segmentize().linestring_idx = np.array([0])
        linestrings.segmentize().interpolate_points.return_value = dummy_points
        nodes = one_linear_object.interpolate_nodes(
            itertools.count(start=2),
            global_dist_calc_points=74.0,
            fixed_nodes=fixed_nodes,
        )

    assert len(nodes) == 0


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
    assert_array_equal(lines.content_type, ContentType.TYPE_V2_1D_LATERAL)
    assert_array_equal(lines.content_pk, [1, 1, 2, 2, 2])
    assert_array_equal(lines.kcu, [2, 2, 1, 1, 1])
    assert_array_equal(lines.cross_id1, [17, 17, 18, 18, 18])
    assert_array_equal(lines.cross_id2, [17, 17, 18, 18, 18])
    assert_array_equal(lines.cross_weight, 1.0)
    assert_array_equal(lines.discharge_coefficient_positive, [0.8, 1.0, 0.5, 1.0, 1.0])
    assert_array_equal(lines.discharge_coefficient_negative, [1.0, 0.5, 1.0, 1.0, 0.8])
    assert_almost_equal(lines.s1d, expected_centers)
    assert_almost_equal(lines.ds1d, expected_sizes)
    assert_almost_equal(shapely.length(lines.line_geometries), expected_sizes)
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
    two_linear_objects.the_geom[0] = shapely.linestrings(geom) if geom else None
    two_linear_objects.set_geometries(connection_nodes)

    expected_geometry = shapely.linestrings(expected)
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


@mock.patch.object(LinearObjects, "has_groundwater_exchange", np.array([False, True]))
def test_apply_has_groundwater_exchange(two_linear_objects):
    nodes = Nodes(
        id=[1, 2, 3, 4],
    )
    lines = Lines(
        id=range(3),
        content_pk=[1, 2, 2],
        content_type=LinearObjects.content_type,
        line=[[2, 1], [1, 3], [3, 4]],
    )
    embedded_nodes = Nodes(id=[])

    two_linear_objects.apply_has_groundwater_exchange(nodes, lines, embedded_nodes)

    assert_array_equal(nodes.has_groundwater_exchange, [1, 0, 1, 1])


@mock.patch.object(LinearObjects, "has_groundwater_exchange", np.array([False, True]))
def test_apply_has_groundwater_exchange_boundary_condition(two_linear_objects):
    nodes = Nodes(
        id=[1, 2, 3, 4],
        calculation_type=[
            CalculationType.CONNECTED,
            CalculationType.CONNECTED,
            CalculationType.CONNECTED,
            CalculationType.BOUNDARY_NODE,
        ],
    )
    lines = Lines(
        id=range(3),
        content_pk=[1, 2, 2],
        content_type=LinearObjects.content_type,
        line=[[2, 1], [1, 3], [3, 4]],
    )
    embedded_nodes = Nodes(id=[])

    two_linear_objects.apply_has_groundwater_exchange(nodes, lines, embedded_nodes)

    assert_array_equal(nodes.has_groundwater_exchange, [1, 0, 1, 0])


@mock.patch.object(LinearObjects, "has_groundwater_exchange", np.array([False, True]))
def test_apply_has_groundwater_exchange_embedded_nodes(two_linear_objects):
    nodes = Nodes(
        id=[1, 2, 3, 4],
        calculation_type=[
            CalculationType.CONNECTED,
            CalculationType.CONNECTED,
            CalculationType.CONNECTED,
            CalculationType.BOUNDARY_NODE,
        ],
    )
    lines = Lines(
        id=range(3),
        content_pk=[1, 2, 2],
        content_type=LinearObjects.content_type,
        line=[[2, 1], [1, 3], [3, 4]],
    )
    embedded_nodes = Nodes(id=[1], embedded_in=[3])

    two_linear_objects.apply_has_groundwater_exchange(nodes, lines, embedded_nodes)

    assert_array_equal(nodes.has_groundwater_exchange, [1, 0, 0, 0])
