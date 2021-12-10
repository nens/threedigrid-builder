from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal
from threedigrid_builder.constants import CrossSectionShape
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import CrossSectionDefinitions
from threedigrid_builder.grid.cross_section_definitions import tabulate_builtin
from threedigrid_builder.grid.cross_section_definitions import (
    tabulate_closed_rectangle,
)
from threedigrid_builder.grid.cross_section_definitions import tabulate_egg
from threedigrid_builder.grid.cross_section_definitions import (
    tabulate_tabulated,
)
from unittest import mock

import numpy as np
import pytest


SHP = CrossSectionShape


@pytest.fixture
def cross_section_definitions():
    return CrossSectionDefinitions(
        id=[1, 3, 9],
        code=["1", "3", "9"],
        shape=[SHP.CIRCLE, SHP.TABULATED_TRAPEZIUM, SHP.TABULATED_RECTANGLE],
        width=["1.22", "3.7 5.0", "5.0 6.0"],
        height=["", "1 2", "1 2"],
    )


def test_convert_multiple(cross_section_definitions):
    table_1 = np.random.random((9, 2))
    table_2 = np.random.random((4, 2))
    with mock.patch.dict(
        "threedigrid_builder.grid.cross_section_definitions.tabulators",
        {
            SHP.CIRCLE: mock.Mock(return_value=(1, 0.1, None)),
            SHP.TABULATED_TRAPEZIUM: mock.Mock(return_value=(5, 15.0, table_1)),
            SHP.TABULATED_RECTANGLE: mock.Mock(return_value=(6, 11.0, table_2)),
        },
    ):
        actual = cross_section_definitions.convert([1, 3, 9])

        assert len(actual) == 3
        assert_array_equal(actual.id, [0, 1, 2])
        assert_array_equal(actual.content_pk, cross_section_definitions.id)
        assert_array_equal(actual.code, cross_section_definitions.code)
        assert_array_equal(actual.shape, [1, 5, 6])  # see mocks
        assert_array_equal(actual.width_1d, [0.1, 15.0, 11.0])  # see mocks

        assert_array_equal(actual.offset, [0, 0, 9])
        assert_array_equal(actual.count, [0, 9, 4])
        assert_array_equal(actual.tables, np.concatenate([table_1, table_2], axis=0))


def test_convert_multiple_filtered(cross_section_definitions):
    with mock.patch.dict(
        "threedigrid_builder.grid.cross_section_definitions.tabulators",
        {
            SHP.CIRCLE: mock.Mock(return_value=(1, 0.1, None)),
        },
    ):
        actual = cross_section_definitions.convert([1])

        assert len(actual) == 1
        assert_array_equal(actual.id, [0])
        assert_array_equal(actual.content_pk, [1])


def test_convert_nonexisting_id(cross_section_definitions):
    with pytest.raises(
        SchematisationError,
        match=r"One or more objects refer to non-existing cross section definitions \[2, 4\].",
    ):
        cross_section_definitions.convert([1, 2, 3, 4])


def test_tabulate_builtin():
    actual = tabulate_builtin("my-shape", "1.52", "ignored")
    assert actual == ("my-shape", 1.52, None)


def test_tabulate_closed_rectangle():
    shape, width_1d, table = tabulate_closed_rectangle("my-shape", "1.52", "5.2")

    assert shape == CrossSectionShape.TABULATED_RECTANGLE
    assert width_1d == 1.52
    assert_almost_equal(table, np.array([[0.0, 1.52], [5.2, 0.0]]))


def test_tabulate_egg():
    shape, width_1d, table = tabulate_egg("my-shape", "1.52", "ignored")

    assert shape == CrossSectionShape.TABULATED_TRAPEZIUM
    assert width_1d == 1.52

    # the expected table is exactly what inpy returns for a width of 1.52
    expected_table = np.array(
        [
            [0.0, 0.0],
            [0.152, 0.584],
            [0.304, 0.816],
            [0.456, 0.99],
            [0.608, 1.128],
            [0.76, 1.242],
            [0.912, 1.336],
            [1.064, 1.41],
            [1.216, 1.468],
            [1.368, 1.506],
            [1.52, 1.52],
            [1.672, 1.504],
            [1.824, 1.442],
            [1.976, 1.31],
            [2.128, 1.038],
            [2.28, 0.0],
        ]
    )
    assert_almost_equal(table, expected_table, decimal=3)


def test_tabulate_tabulated():
    shape, width_1d, table = tabulate_tabulated("my-shape", "1 2 3", "0 1 2")

    assert shape == "my-shape"
    assert width_1d == 3.0  # the max
    assert_almost_equal(table, np.array([[0, 1], [1, 2], [2, 3]], dtype=float))


@pytest.mark.parametrize(
    "shape,width,height",
    [
        (CrossSectionShape.TABULATED_RECTANGLE, "", "1"),
        (CrossSectionShape.TABULATED_RECTANGLE, "1", ""),
        (CrossSectionShape.TABULATED_RECTANGLE, "", ""),
        (CrossSectionShape.TABULATED_RECTANGLE, "1", "1 2"),
        (CrossSectionShape.TABULATED_RECTANGLE, "1 2", "1"),
        (CrossSectionShape.TABULATED_RECTANGLE, "1 1", "2 1"),
    ],
)
def test_tabulate_tabulated_err(shape, width, height):
    with pytest.raises(SchematisationError):
        tabulate_tabulated(shape, width, height)
