import pytest
from numpy.testing import assert_array_equal

from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType, ContentType
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import ExchangeLines


@pytest.mark.parametrize(
    "content_type",
    [
        ContentType.TYPE_V2_CONNECTION_NODES,
        ContentType.TYPE_V2_PIPE,
        ContentType.TYPE_V2_CULVERT,
    ],
)
def test_get_line_mappings_non_channels(content_type):
    nodes = Nodes(
        id=[1],
        content_type=content_type,
        content_pk=[2],
        calculation_type=[CalculationType.CONNECTED],
    )
    exchange_lines = ExchangeLines(id=[1], channel_id=[2])
    line_node_idx, line_exc_id = exchange_lines.get_line_mappings(nodes)

    assert_array_equal(line_node_idx, [0])
    assert_array_equal(line_exc_id, [-9999])


ISO = CalculationType.ISOLATED
C1 = CalculationType.CONNECTED
C2 = CalculationType.DOUBLE_CONNECTED


@pytest.mark.parametrize(
    "node_calc_type,node_content_pk,channel_id,expected_node_idx,expected_exc_id",
    [
        ([], [], [], [], []),
        ([ISO], [1], [], [], []),
        ([C1], [1], [], [0], [-9999]),
        ([C1], [1], [1], [0], [1]),
        ([C1], [2], [1], [0], [-9999]),
        ([C2], [1], [], [0, 0], [-9999, -9999]),
        ([C2], [1], [1], [0, 0], [1, -9999]),
        ([C2], [2], [1], [0, 0], [-9999, -9999]),
        ([C2], [1], [1, 1], [0, 0], [1, 2]),
        ([C2], [1], [1, 2], [0, 0], [1, -9999]),
        ([C1, ISO, C2], [1, 2, 3], [1, 3, 3], [0, 2, 2], [1, 2, 3]),
    ],
)
def test_get_line_mappings_channels(
    node_calc_type, node_content_pk, channel_id, expected_node_idx, expected_exc_id
):
    nodes = Nodes(
        id=range(len(node_content_pk)),
        content_type=ContentType.TYPE_V2_CHANNEL,
        content_pk=node_content_pk,
        calculation_type=node_calc_type,
    )
    exchange_lines = ExchangeLines(
        id=range(1, len(channel_id) + 1), channel_id=channel_id
    )

    line_node_idx, line_exc_id = exchange_lines.get_line_mappings(nodes)

    assert_array_equal(line_node_idx, expected_node_idx)
    assert_array_equal(line_exc_id, expected_exc_id)


def test_get_line_mappings_too_many_exchange_lines():
    exchange_lines = ExchangeLines(id=[1, 2, 3], channel_id=[1, 1, 1])

    with pytest.raises(
        SchematisationError, match=r"Channels \[1\] have too many exchange lines"
    ):
        exchange_lines.get_line_mappings(Nodes(id=[]))
