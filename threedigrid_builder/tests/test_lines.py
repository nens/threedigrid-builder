from numpy.testing import assert_almost_equal
from numpy.testing import assert_equal
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes

import numpy as np
import pygeos
import pytest


@pytest.fixture
def nodes():
    return Nodes(id=[1, 2, 3], coordinates=[(1, 1), (2, 2), (3, 3)])


@pytest.fixture
def lines():
    return Lines(
        id=[0, 1, 2],
        line=[(1, 2), (2, 3), (1, 3)],
        line_geometries=[None, pygeos.linestrings([(5, 5), (6, 6)]), None],
        ds1d=[np.nan, 16.0, 2.0],
    )


def test_set_line_coords(nodes, lines):
    lines.set_line_coords(nodes)

    assert_equal(lines.line_coords, [[1, 1, 2, 2], [2, 2, 3, 3], [1, 1, 3, 3]])


def test_set_line_coords_skips_already_set(nodes, lines):
    lines.line_coords[0] = (9, 9, 10, 10)
    lines.line_coords[1] = (1, 1, 1, np.nan)  # should be recomputed
    lines.set_line_coords(nodes)

    assert_equal(lines.line_coords[:2], [[9, 9, 10, 10], [2, 2, 3, 3]])


def test_fix_line_geometries(lines):
    lines.line_coords = np.array(
        [[1, 1, 2, 2], [np.nan, np.nan, np.nan, np.nan], [1, 1, 3, 3]]
    )

    lines.fix_line_geometries()

    assert pygeos.equals(lines.line_geometries[0], pygeos.linestrings([(1, 1), (2, 2)]))
    assert pygeos.equals(lines.line_geometries[1], pygeos.linestrings([(5, 5), (6, 6)]))
    assert pygeos.equals(lines.line_geometries[2], pygeos.linestrings([(1, 1), (3, 3)]))
    assert_almost_equal(lines.ds1d, [np.sqrt(2), 16.0, 2 * np.sqrt(2)])


def test_set_discharge_coefficients(lines):
    lines.discharge_coefficient_positive[:] = 2.0
    lines.discharge_coefficient_negative[0] = 3.0
    lines.set_discharge_coefficients()

    assert_equal(lines.discharge_coefficient_positive, [2.0, 2.0, 2.0])
    assert_equal(lines.discharge_coefficient_negative, [3.0, 1.0, 1.0])


def test_sort_by_nodes(lines):
    lines.line[:] = [[4, 5], [6, 7], [1, 0]]
    lines.sort_by_nodes([1, 4])

    assert_equal(lines.line, [[1, 0], [6, 7], [4, 5]])
