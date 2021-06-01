from numpy.testing import assert_array_equal
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import ConnectionNodes
from threedigrid_builder.grid import Orifices
from threedigrid_builder.grid import Weirs

import itertools
import numpy as np
import pytest


@pytest.fixture
def connection_nodes():
    # Used to map connection_node_start/end_id to an index (sequence id)
    return ConnectionNodes(id=np.array([21, 25, 33, 42]))


@pytest.fixture(scope="module", params=[Weirs, Orifices])
def one_weir_orifice(request):
    return request.param(
        id=[1],
        connection_node_start_id=[21],
        connection_node_end_id=[42],
        crest_level=[2.3],
        crest_type=[4],
    )


@pytest.fixture(scope="module", params=[Weirs, Orifices])
def two_weir_orifices(request):
    return request.param(
        id=[1, 2],
        connection_node_start_id=[21, 25],
        connection_node_end_id=[42, 33],
        crest_level=[2.3, 4.5],
        crest_type=[4, 3],
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
    assert_array_equal(lines.line, [(100, 101), (103, 102)])
    assert_array_equal(lines.content_type, expected_content_type)
    assert_array_equal(lines.content_pk, [1, 2])
    assert_array_equal(lines.kcu, [4, 3])
    assert_array_equal(lines.dpumax, [2.3, 4.5])
