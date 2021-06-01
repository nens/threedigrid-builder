from numpy.testing import assert_array_equal
from threedigrid_builder.base import Nodes
from threedigrid_builder.grid import Channels
from threedigrid_builder.grid import CrossSectionLocations

import numpy as np
import pygeos
import pytest


@pytest.fixture
def channels():
    return Channels(
        the_geom=pygeos.linestrings([[(0, 0), (6, 0), (6, 6)]]),
        dist_calc_points=np.array([5.0]),
        id=np.array([1]),
        code=np.array(["one"]),
        connection_node_start_id=np.array([21]),
        connection_node_end_id=np.array([42]),
        calculation_type=np.array([2]),
    )


def test_1d2d_properties(channels):
    nodes = Nodes(
        id=[0, 2, 5, 7],
        dmax=[1.0, 3.0, 2.0, 4.0],
        cross1=[2, 5, -9999, 10],
        cross2=[5, 10, -9999, 13],
        cross_weight=[0.2, 0.5, np.nan, 0.8],
    )
    locations = CrossSectionLocations(
        id=[2, 5, 10, 13],
        bank_level=[1.0, 2.0, 3.0, 4.0],
    )
    node_idx = [0, 1, 3]

    is_sewerage, dpumax = channels.get_1d2d_properties(nodes, node_idx, locations)

    # channels are no sewerage
    assert_array_equal(is_sewerage, False)

    # for the manholes, conn_node.drain_level is copied, otherwise, dmax is taken
    assert_array_equal(dpumax, [1.8, 2.5, 3.2])
