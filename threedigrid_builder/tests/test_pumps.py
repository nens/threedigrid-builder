import numpy as np
import pytest
from numpy.testing import assert_array_equal

from threedigrid_builder.base import Nodes, Pumps
from threedigrid_builder.constants import ContentType
from threedigrid_builder.exceptions import SchematisationError


@pytest.fixture
def nodes():
    CN = ContentType.TYPE_V2_CONNECTION_NODES
    CH = ContentType.TYPE_V2_CHANNEL
    return Nodes(
        id=[0, 1, 3, 4],
        content_type=[CH, CN, CN, CH],
        content_pk=[21, 21, 33, 42],
        coordinates=[(0, 21), (1, 25), (2, 33), (3, 42)],
        dmax=[1.0, 2.0, 3.0, 4.0],
    )


@pytest.fixture
def pumps():
    return Pumps(
        id=[2, 35],
        line=[(1, 3), (1, -9999)],
        connection_node_start_id=[21, 21],
        connection_node_end_id=[33, -9999],
    )


def test_renumber(pumps):
    pumps.renumber()
    assert_array_equal(pumps.id, [0, 1])
    assert_array_equal(pumps.content_pk, [2, 35])


def test_set_lines(pumps, nodes):
    pumps.line[:] = -9999
    pumps.set_lines(nodes)
    assert_array_equal(pumps.line, [(1, 3), (1, -9999)])


@pytest.mark.parametrize(
    "start,end",
    [
        ([-9999, 21], [33, 33]),
        ([42, 21], [33, 33]),
        ([21, 21], [33, 42]),
    ],
)
def test_set_lines_error(start, end, pumps, nodes):
    pumps.connection_node_start_id[:] = start
    pumps.connection_node_end_id[:] = end
    with pytest.raises(SchematisationError):
        pumps.set_lines(nodes)


def test_set_node_data(pumps, nodes):
    pumps.set_node_data(nodes)
    assert_array_equal(pumps.line_coords, [(1, 25, 2, 33), (1, 25, np.nan, np.nan)])
    assert_array_equal(pumps.coordinates, [(1.5, 29.0), (1.0, 25.0)])
    assert_array_equal(pumps.bottom_level, [2.0, 2.0])
