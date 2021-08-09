from numpy.testing import assert_array_equal
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import ContentType
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
        id=[0, 2, 5],
        s1d=[3, np.nan, 6],
        content_type=[
            ContentType.TYPE_V2_CHANNEL,
            ContentType.TYPE_V2_CONNECTION_NODES,
            ContentType.TYPE_V2_CHANNEL,
        ],
        content_pk=[1, 3, 1],
    )
    locations = CrossSectionLocations(
        id=[2, 5],
        the_geom=pygeos.points([(0, 0), [6, 6]]),
        bank_level=[1.0, 13.0],
        channel_id=[1, 1],
    )
    node_idx = [0, 2]

    is_closed, dpumax = channels.get_1d2d_properties(nodes, node_idx, locations)

    # channels are no sewerage
    assert_array_equal(is_closed, False)

    # bank levels are interpolated
    assert_array_equal(dpumax, [4.0, 7.0])
