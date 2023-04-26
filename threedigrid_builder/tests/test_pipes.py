import itertools
from unittest import mock

import numpy as np
import pytest
import shapely
from numpy.testing import assert_almost_equal, assert_array_equal

from threedigrid_builder.base import Nodes
from threedigrid_builder.grid import ConnectionNodes, Pipes


@pytest.fixture
def connection_nodes():
    # Used to map connection_node_start/end_id to an index (sequence id)
    return ConnectionNodes(
        id=[21, 25, 33, 42],
        the_geom=shapely.points([(0, 21), (1, 25), (2, 33), (3, 42)]),
        drain_level=[0.0, 20.0, 30.0, 10.0],
    )


@pytest.fixture
def pipes():
    return Pipes(
        id=[1, 2],
        dist_calc_points=[5.0, np.nan],
        code=["one", "two"],
        connection_node_start_id=[21, 21],
        connection_node_end_id=[25, 42],
        calculation_type=[2, 1],
        invert_level_start_point=[3.0, 5.0],
        invert_level_end_point=[4.0, 6.0],
    )


@pytest.mark.parametrize(
    "pipe_ids,ds,expected",
    [
        ([1], [5.0], [3.5]),
        ([2], [0.5], [5.5]),
        ([1, 2], [5.0, 0.5], [3.5, 5.5]),
        ([2, 1], [0.5, 5.0], [5.5, 3.5]),
        ([1, 1, 2, 2], [5.0, 7.5, 0.5, 0.25], [3.5, 3.75, 5.5, 5.25]),
    ],
)
def test_compute_bottom_level(pipe_ids, ds, pipes, expected):
    # set geometries with lengths 10 and 1 (resp. id 1 and 2)
    # invert levels are [3, 4] for id=1 and [5, 6] for id=2
    pipes.the_geom = shapely.linestrings([[(0, 0), (0, 10)], [(2, 2), (3, 2)]])

    actual = pipes.compute_bottom_level(pipe_ids, ds)

    assert_almost_equal(actual, expected)


def test_compute_bottom_level_raises_no_geom(pipes):
    with pytest.raises(ValueError, match=".*Call set_geometries first.*"):
        pipes.compute_bottom_level([1], [5.0])


@pytest.mark.parametrize(
    "pipe_ids,ds",
    [
        ([1], [10.1]),
        ([2], [-1e-7]),
    ],
)
def test_compute_bottom_level_raises_out_of_bounds(pipe_ids, ds, pipes):
    pipes.the_geom = shapely.linestrings([[(0, 0), (0, 10)], [(2, 2), (3, 2)]])
    with pytest.raises(ValueError, match=".*outside of the linear object bounds.*"):
        pipes.compute_bottom_level(pipe_ids, ds)


def test_get_1d2d_exchange_levels(pipes):
    content_pk = np.array([1, 1, 2])
    s1d = np.array([12.0, 13.0, 15.0])
    connection_nodes = mock.Mock()

    with mock.patch.object(pipes, "compute_drain_level") as compute_drain_level:
        compute_drain_level.return_value = np.array([1.8, 2.5, 3.2])
        actual = pipes.get_1d2d_exchange_levels(
            content_pk=content_pk, s1d=s1d, connection_nodes=connection_nodes
        )

        _, kwargs = compute_drain_level.call_args
        assert_array_equal(kwargs["ids"], content_pk)
        assert_array_equal(kwargs["s"], s1d)
        assert kwargs["connection_nodes"] is connection_nodes

    # interpolation between manhole drain levels is further tested elsewhere
    assert_array_equal(actual, compute_drain_level.return_value)


def test_is_closed(pipes):
    actual = pipes.is_closed(np.array([1, 3]))

    assert_array_equal(actual, [True, True])


@pytest.mark.parametrize(
    "thickness,hc_out,hc_in,expected",
    [
        (0.1, 3.0, 2.0, True),
        (0.0, 3.0, 2.0, False),
        (np.nan, 3.0, 2.0, False),
        (0.1, np.nan, 2.0, False),
        (0.1, 3.0, np.nan, False),
    ],
)
def test_has_groundwater_exchange(thickness, hc_out, hc_in, expected):
    pipes = Pipes(
        id=[1],
        exchange_thickness=thickness,
        hydraulic_conductivity_out=hc_out,
        hydraulic_conductivity_in=hc_in,
    )

    actual = pipes.has_groundwater_exchange

    assert len(actual) == 1
    assert actual[0] == expected


@pytest.mark.parametrize(
    "thickness,hc_out,hc_in,res_out,res_in",
    [
        (0.1, 3.0, 2.0, 30.0, 20.0),
        (np.nan, 3.0, 2.0, np.nan, np.nan),
        (0.1, np.nan, 2.0, np.nan, 20.0),
        (0.1, 3.0, np.nan, 30.0, np.nan),
    ],
)
def test_set_hydraulic_resistance(thickness, hc_out, hc_in, res_in, res_out):
    pipes = Pipes(
        the_geom=shapely.linestrings([[(0, 0), (6, 0), (6, 6)]]),
        dist_calc_points=np.array([5.0]),
        id=np.array([1]),
        code=np.array(["one"]),
        connection_node_start_id=np.array([21]),
        connection_node_end_id=np.array([42]),
        calculation_type=np.array([2]),
    )
    pipes.exchange_thickness = np.array([thickness])
    pipes.hydraulic_conductivity_out = np.array([hc_out])
    pipes.hydraulic_conductivity_in = np.array([hc_in])

    cn = ConnectionNodes(id=[21, 42], the_geom=shapely.points([(0, 21), (1, 25)]))
    nodes = Nodes(id=[10], content_pk=[1], s1d=[0.5])
    lines = pipes.get_lines(cn, nodes, itertools.count(start=0))

    assert_array_equal(lines.hydraulic_resistance_in, [res_in, res_in])
    assert_array_equal(lines.hydraulic_resistance_out, [res_out, res_out])
