import itertools
from unittest import mock

import numpy as np
import pytest
import shapely
from numpy.testing import assert_almost_equal, assert_array_equal

from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import ConnectionNodes, Culverts, Orifices, Weirs


@pytest.fixture
def connection_nodes():
    # Used to map connection_node_start/end_id to an index (sequence id)
    return ConnectionNodes(
        id=np.array([21, 25, 33, 42]),
        the_geom=shapely.points([(0, 21), (1, 25), (2, 33), (3, 42)]),
    )


@pytest.fixture
def culverts():
    return Culverts(
        id=[1, 2],
        dist_calc_points=[5.0, np.nan],
        code=["one", "two"],
        connection_node_start_id=[21, 21],
        connection_node_end_id=[25, 42],
        calculation_type=[2, 1],
        invert_level_start_point=[3.0, 5.0],
        invert_level_end_point=[4.0, 6.0],
        the_geom=[
            shapely.linestrings([(0, 21), (0.5, 22), (1, 25)]),
            shapely.linestrings([(0, 21), (3, 42)]),
        ],
    )


@pytest.fixture(scope="module", params=[Weirs, Orifices])
def two_weir_orifices(request):
    return request.param(
        id=[1, 2],
        connection_node_start_id=[21, 25],
        connection_node_end_id=[42, 33],
        crest_level=[2.3, 4.5],
        crest_type=[4, 3],
        friction_type=[1, 1],
        friction_value=[31, 41],
        cross_section_definition_id=[17, 18],
    )


def test_get_lines(connection_nodes, two_weir_orifices):
    lines = two_weir_orifices.get_lines(
        connection_nodes,
        itertools.count(start=0),
        connection_node_offset=100,
    )

    if two_weir_orifices.__class__ is Weirs:
        expected_content_type = ContentType.TYPE_V2_WEIR
    elif two_weir_orifices.__class__ is Orifices:
        expected_content_type = ContentType.TYPE_V2_ORIFICE

    assert_array_equal(lines.id, range(2))
    assert_array_equal(lines.line, [(100, 103), (101, 102)])
    assert_array_equal(lines.content_type, expected_content_type)
    assert_array_equal(lines.content_pk, [1, 2])
    assert_array_equal(lines.kcu, [4, 3])
    assert_array_equal(lines.cross_id1, [17, 18])
    assert_array_equal(lines.cross_id2, [17, 18])
    assert_array_equal(lines.cross_weight, [1.0, 1.0])
    assert_array_equal(lines.dpumax, [2.3, 4.5])
    assert_array_equal(lines.invert_level_start_point, [2.3, 4.5])
    assert_array_equal(lines.invert_level_end_point, [2.3, 4.5])
    assert_array_equal(lines.frict_type1, [1, 1])
    assert_array_equal(lines.frict_type2, [1, 1])
    assert_array_equal(lines.frict_value1, [31, 41])
    assert_array_equal(lines.frict_value2, [31, 41])


@pytest.mark.parametrize(
    "culvert_ids,ds,expected",
    [
        ([1], [5.0], [3.5]),
        ([2], [0.5], [5.5]),
        ([1, 2], [5.0, 0.5], [3.5, 5.5]),
        ([2, 1], [0.5, 5.0], [5.5, 3.5]),
        ([1, 1, 2, 2], [5.0, 7.5, 0.5, 0.25], [3.5, 3.75, 5.5, 5.25]),
    ],
)
def test_culverts_compute_bottom_level(culvert_ids, ds, culverts, expected):
    # set geometries with lengths 10 and 1 (resp. id 1 and 2)
    # invert levels are [3, 4] for id=1 and [5, 6] for id=2
    culverts.the_geom = shapely.linestrings([[(0, 0), (0, 10)], [(2, 2), (3, 2)]])

    actual = culverts.compute_bottom_level(culvert_ids, ds)

    assert_almost_equal(actual, expected)


def test_culverts_compute_bottom_level_raises_no_geom(culverts):
    culverts.the_geom[:] = None
    with pytest.raises(ValueError, match=".*Call set_geometries first.*"):
        culverts.compute_bottom_level([1], [0.2])


@pytest.mark.parametrize(
    "culvert_ids,ds",
    [
        ([1], [10.1]),
        ([2], [-1e-7]),
    ],
)
def test_culverts_compute_bottom_level_raises_out_of_bounds(culvert_ids, ds, culverts):
    culverts.the_geom = shapely.linestrings([[(0, 0), (0, 10)], [(2, 2), (3, 2)]])
    with pytest.raises(ValueError, match=".*outside of the linear object bounds.*"):
        culverts.compute_bottom_level(culvert_ids, ds)


def test_culverts_get_1d2d_exchange_levels(culverts):
    content_pk = np.array([1, 1, 2])
    s1d = np.array([12.0, 13.0, 15.0])
    connection_nodes = mock.Mock()

    with mock.patch.object(culverts, "compute_drain_level") as compute_drain_level:
        compute_drain_level.return_value = np.array([1.8, 2.5, 3.2])
        actual = culverts.get_1d2d_exchange_levels(
            content_pk=content_pk, s1d=s1d, connection_nodes=connection_nodes
        )

        _, kwargs = compute_drain_level.call_args
        assert_array_equal(kwargs["ids"], content_pk)
        assert_array_equal(kwargs["s"], s1d)
        assert kwargs["connection_nodes"] is connection_nodes

    # interpolation between manhole drain levels is further tested elsewhere
    assert_array_equal(actual, compute_drain_level.return_value)


def test_is_closed(culverts):
    actual = culverts.is_closed(np.array([1, 3]))

    assert_array_equal(actual, [True, True])
