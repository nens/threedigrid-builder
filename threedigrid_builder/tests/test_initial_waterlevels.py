from unittest import mock

import pytest
import shapely
from numpy import nan
from numpy.testing import assert_almost_equal

from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import Channels, ConnectionNodes
from threedigrid_builder.grid.initial_waterlevels import compute_initial_waterlevels


@pytest.fixture
def channels():
    return Channels(
        id=[51, 52, 53, 54],
        connection_node_start_id=[1, 1, 1, 1],
        connection_node_end_id=[2, 3, 4, 5],
        the_geom=[
            shapely.linestrings([(0, 0), (3, 0)]),  # 1 segment of size 3
            shapely.linestrings([(55, 3), (60, 3)]),  # 1 segment of size 5
            shapely.linestrings([(3, 1), (36, 1)]),  # 3 segments of size 11
            shapely.linestrings([(40, 2), (54, 2)]),  # 2 segments of size 7
        ],
    )


@pytest.fixture
def nodes():
    CH = ContentType.TYPE_V2_CHANNEL
    CN = ContentType.TYPE_V2_CONNECTION_NODES
    return Nodes(
        id=range(8),
        content_pk=[1, 2, 3, 53, 53, 4, 54, 5],
        content_type=[CN, CN, CN, CH, CH, CN, CH, CN],
        s1d=[nan, nan, nan, 11.0, 22.0, nan, 7.0, nan],
        dmax=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8],
    )


@pytest.fixture
def connection_nodes():
    return ConnectionNodes(id=[1, 2, 3, 4, 5])


@pytest.mark.parametrize(
    "initial_waterlevels, expected",
    [
        ([0.0, 0.0, 0.0, 0.0, 0.0], [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]),
        ([1.0, 1.0, 1.0, 1.0, 1.0], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
        ([1.0, nan, nan, nan, nan], [1.0, 0.2, 0.3, 0.867, 0.733, 0.6, 0.9, 0.8]),
        ([nan, nan, nan, nan, 1.0], [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0]),
        ([nan, nan, nan, nan, 2.0], [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.05, 2.0]),
    ],
)
def test_compute_initial_waterlevels(
    nodes, channels, connection_nodes, initial_waterlevels, expected
):
    connection_nodes.initial_waterlevel[:] = initial_waterlevels
    compute_initial_waterlevels(
        nodes,
        connection_nodes=connection_nodes,
        channels=channels,
        pipes=mock.Mock(),
        culverts=mock.Mock(),
    )

    assert_almost_equal(nodes.initial_waterlevel, expected, decimal=3)
