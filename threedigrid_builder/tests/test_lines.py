import numpy as np
import pytest
import shapely
from numpy.testing import assert_almost_equal, assert_equal
from shapely.testing import assert_geometries_equal

from threedigrid_builder.base import Lines, Nodes
from threedigrid_builder.constants import LineType


@pytest.fixture
def nodes():
    return Nodes(id=[1, 2, 3], coordinates=[(1, 1), (2, 2), (3, 3)])


@pytest.fixture
def lines():
    return Lines(
        id=[0, 1, 2],
        line=[(1, 2), (2, 3), (1, 3)],
        line_geometries=[None, shapely.linestrings([(5, 5), (6, 6)]), None],
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

    expected = shapely.linestrings(
        [[(1, 1), (2, 2)], [(5, 5), (6, 6)], [(1, 1), (3, 3)]]
    )
    assert_geometries_equal(lines.line_geometries, expected)


def test_fix_ds1d(lines):
    lines.line_geometries = shapely.linestrings(
        [[(1, 1), (2, 2)], [(5, 5), (6, 6)], [(1, 1), (3, 3)]]
    )
    lines.fix_ds1d()

    assert_almost_equal(lines.ds1d, [np.sqrt(2), 16.0, 2.0])


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


U = LineType.LINE_2D_U
V = LineType.LINE_2D_V
U_OBS = LineType.LINE_2D_OBSTACLE_U
V_OBS = LineType.LINE_2D_OBSTACLE_V
L1D2D = LineType.LINE_1D2D_SINGLE_CONNECTED_CLOSED


@pytest.mark.parametrize(
    "kcu,crest_levels,expected_kcu,expected_flod",
    [
        ([U, V, V], [0.1, 0.3], [U_OBS, V_OBS, V], [0.1, 0.3, np.nan]),
        ([U, U, V], [], [U, U, V], [np.nan, np.nan, np.nan]),
    ],
)
def test_set_2d_crest_levels(
    lines: Lines, kcu, crest_levels, expected_kcu, expected_flod
):
    lines.kcu[:] = kcu
    lines.set_2d_crest_levels(
        np.array(crest_levels), where=np.arange(len(crest_levels))
    )

    assert_equal(lines.kcu, expected_kcu)
    assert_equal(lines.flod, expected_flod)
    assert_equal(lines.flou, expected_flod)


def test_line_get_velocity_points(lines: Lines):
    lines.ds1d_half[1] = np.sqrt(2) * 0.5

    actual = lines.get_velocity_points([1])
    assert len(actual) == 1
    assert shapely.to_wkt(actual[0]) == "POINT (5.5 5.5)"
