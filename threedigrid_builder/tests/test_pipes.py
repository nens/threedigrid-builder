from numpy.testing import assert_almost_equal
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import Pipes

import numpy as np
import pygeos
import pytest


@pytest.fixture
def connection_nodes():
    # Used to map connection_node_start/end_id to an index (sequence id)
    return ConnectionNodes(
        id=[21, 25, 33, 42],
        the_geom=pygeos.points([(0, 21), (1, 25), (2, 33), (3, 42)]),
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
    pipes.the_geom = pygeos.linestrings([[(0, 0), (0, 10)], [(2, 2), (3, 2)]])

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
    pipes.the_geom = pygeos.linestrings([[(0, 0), (0, 10)], [(2, 2), (3, 2)]])
    with pytest.raises(ValueError, match=".*outside of the linear object bounds.*"):
        pipes.compute_bottom_level(pipe_ids, ds)
