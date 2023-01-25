import numpy as np
import pytest
from numpy.testing import assert_array_equal

from threedigrid_builder.grid import ExchangeLines


def test_channel_mapping():
    exchange_lines = ExchangeLines(id=range(1, 6), channel_id=[15, 2, 3, 15, 3])
    actual = exchange_lines.channel_mapping

    assert_array_equal(actual.id, [2, 3, 15])
    assert_array_equal(actual.exchange_line_id, [2, 3, 1])
    assert_array_equal(actual.secondary_exchange_line_id, [-9999, 5, 4])


@pytest.mark.parametrize(
    "channel_id,expected,is_primary",
    [
        (np.array([]), [], True),
        (np.array([1]), [-9999], True),
        (np.array([2]), [2], True),
        (np.array([2]), [-9999], False),
        (np.array([3]), [3], True),
        (np.array([3]), [5], False),
    ],
)
def test_get_for_channel_id(channel_id, expected, is_primary):
    exchange_lines = ExchangeLines(id=range(1, 6), channel_id=[15, 2, 3, 15, 3])
    actual = exchange_lines.get_for_channel_id(channel_id, is_primary)

    assert_array_equal(actual, expected)
