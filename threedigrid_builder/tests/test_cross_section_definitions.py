from unittest import mock

import numpy as np
import pytest
from numpy.testing import assert_almost_equal, assert_array_equal
from threedi_schema import constants

from threedigrid_builder.constants import CrossSectionShape
from threedigrid_builder.exceptions import SchematisationError
from threedigrid_builder.grid import CrossSectionDefinitions, CrossSections
from threedigrid_builder.grid.cross_section_definitions import (
    _parse_tabulated,
    set_friction_vegetation_values,
    tabulate_builtin,
    tabulate_closed_rectangle,
    tabulate_egg,
    tabulate_inverted_egg,
    tabulate_tabulated,
    tabulate_yz,
)

SHP = CrossSectionShape


@pytest.fixture
def cross_section_definitions():
    return CrossSectionDefinitions(
        id=[1, 3, 9],
        code=["1", "3", "9"],
        shape=[SHP.CIRCLE, SHP.TABULATED_TRAPEZIUM, SHP.TABULATED_RECTANGLE],
        width=[1.22, None, None],
        height=[None, None, None],
        cross_section_table=[None, "1,3.7\n2,5.0", "1,5\n2,6"],
    )


def test_convert_multiple(cross_section_definitions):
    table_1 = np.random.random((9, 2))
    table_2 = np.random.random((4, 2))
    with mock.patch.dict(
        "threedigrid_builder.grid.cross_section_definitions.tabulators",
        {
            SHP.CIRCLE: mock.Mock(return_value=(1, 0.1, None, None, None)),
            SHP.TABULATED_TRAPEZIUM: mock.Mock(
                return_value=(5, 15.0, 2.0, table_1, None)
            ),
            SHP.TABULATED_RECTANGLE: mock.Mock(
                return_value=(6, 11.0, 2.0, table_2, None)
            ),
        },
    ):
        actual = cross_section_definitions.convert([1, 3, 9])

        assert len(actual) == 3
        assert_array_equal(actual.id, [0, 1, 2])
        assert_array_equal(actual.content_pk, cross_section_definitions.id)
        assert_array_equal(actual.code, cross_section_definitions.code)
        assert_array_equal(actual.shape, [1, 5, 6])  # see mocks
        assert_array_equal(actual.width_1d, [0.1, 15.0, 11.0])  # see mocks
        assert_array_equal(actual.height_1d, [np.nan, 2.0, 2.0])
        assert_array_equal(actual.offset, [0, 0, 9])
        assert_array_equal(actual.count, [0, 9, 4])
        assert_array_equal(actual.tables, np.concatenate([table_1, table_2], axis=0))


def test_convert_multiple_filtered(cross_section_definitions):
    with mock.patch.dict(
        "threedigrid_builder.grid.cross_section_definitions.tabulators",
        {
            SHP.CIRCLE: mock.Mock(return_value=(1, 0.1, 0.1, None, None)),
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


def test_convert_empty(cross_section_definitions):
    result = cross_section_definitions.convert([])
    assert isinstance(result, CrossSections)
    assert len(result) == 0


def test_tabulate_builtin():
    actual = tabulate_builtin("my-shape", 1.52, 1.33)
    assert actual == ("my-shape", 1.52, None, None, None)


def test_tabulate_closed_rectangle():
    shape, width_1d, height_1d, table, yz = tabulate_closed_rectangle(
        "my-shape", 1.52, 5.2
    )

    assert shape == CrossSectionShape.TABULATED_RECTANGLE
    assert width_1d == 1.52
    assert height_1d == 5.2
    assert_almost_equal(table, np.array([[0.0, 1.52], [5.2, 0.0]]))
    assert yz is None


def test_tabulate_egg():
    shape, width_1d, height_1d, table, yz = tabulate_egg("my-shape", 1.52, 0)

    assert shape == CrossSectionShape.TABULATED_TRAPEZIUM
    assert width_1d == 1.52
    assert height_1d == 1.52 * 1.5
    assert yz is None

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
    width, height = _parse_tabulated(
        "0,1\n1,2\n2,3", CrossSectionShape.TABULATED_RECTANGLE
    )
    shape, width_1d, height_1d, table, yz = tabulate_tabulated(
        "my-shape", width, height
    )
    assert shape == "my-shape"
    assert width_1d == 3.0  # the max
    assert height_1d == 2.0
    assert_almost_equal(table, np.array([[0, 1], [1, 2], [2, 3]], dtype=float))
    assert yz is None


def test_tabulate_tabulated_err():
    with pytest.raises(
        SchematisationError, match=r".*of tabulated type must have increasing heights.*"
    ):
        tabulate_tabulated(CrossSectionShape.TABULATED_RECTANGLE, [1, 1], [2, 1])


def test_tabulate_inverted_egg():
    shape, width_1d, height_1d, table, yz = tabulate_inverted_egg(
        "my-shape", 1.52, "ignored"
    )

    assert shape == CrossSectionShape.TABULATED_TRAPEZIUM
    assert width_1d == 1.52
    assert height_1d == 1.52 * 1.5
    assert yz is None

    # the expected table is exactly what inpy returns for a width of 1.52
    expected_table = np.array(
        [
            [0.0, 0.0],
            [0.152, 1.038],
            [0.304, 1.31],
            [0.456, 1.442],
            [0.608, 1.504],
            [0.76, 1.52],
            [0.912, 1.506],
            [1.064, 1.468],
            [1.216, 1.41],
            [1.368, 1.336],
            [1.52, 1.242],
            [1.672, 1.128],
            [1.824, 0.99],
            [1.976, 0.816],
            [2.128, 0.584],
            [2.28, 0.0],
        ]
    )
    assert_almost_equal(table, expected_table, decimal=3)


@pytest.mark.parametrize(
    "shape, expected_width, expected_height",
    [
        (CrossSectionShape.TABULATED_YZ, [1, 3, 5], [2, 4, 6]),
        (CrossSectionShape.TABULATED_RECTANGLE, [2, 4, 6], [1, 3, 5]),
    ],
)
def test_parse_tabulated(shape, expected_width, expected_height):
    width, height = _parse_tabulated("1,2\n3,4\n5,6", shape)
    np.testing.assert_array_equal(width, expected_width)
    np.testing.assert_array_equal(height, expected_height)


@pytest.mark.parametrize("cross_section_table", ["", "1,2\n3", "1,2\n3,"])
def test_parse_tabulated_invalid(cross_section_table):
    with pytest.raises(SchematisationError):
        _parse_tabulated(cross_section_table, CrossSectionShape.TABULATED_RECTANGLE)


@pytest.mark.parametrize(
    "cross_section_table,friction_values,vegetation_stem_densities,vegetation_stem_diameters,vegetation_heights,vegetation_drag_coefficients,exp_width,exp_height,exp_table,exp_yz",
    [
        (
            "0,0.5\n0.5,0\n1,0\n1.5,0.5",
            None,
            None,
            None,
            None,
            None,
            1.5,
            0.5,
            [[0, 0.5], [0.5, 1.5]],
            [
                [0, 0.5, 0, 0, 0],
                [0.5, 0, 0, 0, 0],
                [1, 0, 0, 0, 0],
                [1.5, 0.5, 0, 0, 0],
            ],
        ),
        (
            "0,0.5\n0.5,0\n1,0\n1.5,0.25",
            "1 1 1",
            "0.5 1 1",
            "1 1 1",
            "3 1 1",
            "2 1 1",
            1.5,
            0.5,
            [[0, 0.5], [0.25, 1.25], [0.5, 1.5]],
            [
                [0, 0.5, 1, 1, 3],
                [0.5, 0, 1, 1, 1],
                [1, 0, 1, 1, 1],
                [1.5, 0.25, 0, 0, 0],
            ],
        ),
        (
            "0,1\n1,0\n2,0.5\n3,0.5\n4,0\n5,1",
            "1 1 1 1 1",
            "1 1 1 1 1",
            "1 1 1 1 1",
            "1 1 1 1 1",
            "1 1 1 1 1",
            5,
            1,
            [[0, 0], [0.5, 3], [0.5, 4], [1, 5]],
            [
                [0, 1, 1, 1, 1],
                [1, 0, 1, 1, 1],
                [2, 0.5, 1, 1, 1],
                [3, 0.5, 1, 1, 1],
                [4, 0, 1, 1, 1],
                [5, 1, 0, 0, 0],
            ],
        ),
        (
            "0,0.5\n1,0\n2,0.5\n2,1.5\n0,1.5\n0,0.5",
            "1 1 1 1",
            "1 1 1 1",
            "1 1 1 1",
            "1 1 1 1",
            "1 1 1 1",
            2.0,
            1.5,
            [[0, 0], [0.5, 2.0], [1.5, 2.0], [1.5, 0.0]],
            None,
        ),
        (
            "0,0.5\n0.5,0\n0.75,0\n1.0,0\n1.5,0.5",
            "1 1 1 1",
            "1 0.1 1 1",
            "2 1 1 1",
            "1 1 1 1",
            "1 1 1 1",
            1.5,
            0.5,
            [[0, 0.5], [0.5, 1.5]],
            [
                [0, 0.5, 1, 2, 1],
                [0.5, 0, 1, 0.1, 1],
                [0.75, 0, 1, 1, 1],
                [1, 0, 1, 1, 1],
                [1.5, 0.5, 0, 0, 0],
            ],
        ),
        (
            "0,0\n1,1\n0,1\n1,0\n0,0",
            "1 1 1 1",
            "1 1 1 1",
            "1 1 1 1",
            "1 1 1 1",
            "1 1 1 1",
            1,
            1,
            [[0, 1], [0.5, 0], [1, 1], [1, 0]],
            None,
        ),
        (  # Open profile left side higher than right side
            "0,1\n1,0\n2,0\n3,0\n4,2",
            None,
            None,
            None,
            None,
            None,
            4.0,
            2.0,
            np.column_stack(
                (
                    [0.0, 1.0, 2.0],
                    [2.0, 3.5, 4.0],
                )
            ),
            np.column_stack(
                (
                    [0.0, 1.0, 2.0, 3.0, 4.0],
                    [1.0, 0.0, 0.0, 0.0, 2.0],
                    np.zeros(5, dtype=float),
                    np.zeros(5, dtype=float),
                    np.zeros(5, dtype=float),
                )
            ),
        ),
        (  # Open profile right side higher than left side
            "1,3\n2,1\n3,0\n4,1\n5,2",
            None,
            None,
            None,
            None,
            None,
            4.0,
            3.0,
            np.column_stack(
                (
                    [0.0, 1.0, 2.0, 3.0],
                    [0.0, 2.0, 3.5, 4.0],
                )
            ),
            np.column_stack(
                (
                    [1.0, 2.0, 3.0, 4.0, 5.0],
                    [3.0, 1.0, 0.0, 1.0, 2.0],
                    np.zeros(5, dtype=float),
                    np.zeros(5, dtype=float),
                    np.zeros(5, dtype=float),
                )
            ),
        ),
        (  # Open profile same height left and right
            "1,3\n2,1\n3,0\n4,1\n5,3",
            None,
            None,
            None,
            None,
            None,
            4.0,
            3.0,
            np.column_stack(
                (
                    [0.0, 1.0, 3.0],
                    [0.0, 2.0, 4.0],
                )
            ),
            np.column_stack(
                (
                    [1.0, 2.0, 3.0, 4.0, 5.0],
                    [3.0, 1.0, 0.0, 1.0, 3.0],
                    np.zeros(5, dtype=float),
                    np.zeros(5, dtype=float),
                    np.zeros(5, dtype=float),
                )
            ),
        ),
        (  # Open profile, left side rises then falls
            "0,5.2\n5,5.26\n10,5.25\n15,5.2\n20,5.1\n25,5.1\n30,0\n35,5\n40,5",
            None,
            None,
            None,
            None,
            None,
            40.0,
            5.26,
            np.column_stack(
                (
                    [0.0, 5.0, 5.0, 5.1, 5.1, 5.2, 5.2, 5.25, 5.26],
                    [0.0, 9.902, 14.902, 15.0, 20.005, 24.995, 25.008, 34.167, 40.0],
                )
            ),
            np.column_stack(
                (
                    [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0],
                    [5.2, 5.26, 5.25, 5.2, 5.1, 5.1, 0.0, 5.0, 5.0],
                    np.zeros(9, dtype=float),
                    np.zeros(9, dtype=float),
                    np.zeros(9, dtype=float),
                )
            ),
        ),
    ],
)
def test_tabulate_yz(
    cross_section_table,
    friction_values,
    vegetation_stem_densities,
    vegetation_stem_diameters,
    vegetation_heights,
    vegetation_drag_coefficients,
    exp_width,
    exp_height,
    exp_table,
    exp_yz,
):
    width, height = _parse_tabulated(
        cross_section_table, CrossSectionShape.TABULATED_YZ
    )
    shape, width_1d, height_1d, table, yz = tabulate_yz(
        "my-shape",
        width,
        height,
    )

    assert width_1d == exp_width
    assert height_1d == exp_height
    assert_almost_equal(table, np.array(exp_table, dtype=float))
    if yz is not None:
        assert shape == CrossSectionShape.TABULATED_YZ
    else:
        assert shape == CrossSectionShape.TABULATED_TRAPEZIUM


class TestSetFrictionVegetationValues:
    def test_set_friction_values(self):
        friction_values = "0.5, 0.6"
        yz = np.zeros((3, 5), dtype=float)
        result = set_friction_vegetation_values(
            yz, friction_values, None, None, None, None, None
        )
        expected = np.zeros_like(yz)
        expected[:, 2] = [0.5, 0.6, 0]
        np.testing.assert_array_equal(result, expected)

    def test_set_with_cross_section_vegetation_table(self):
        vegetation_table = "1,2,3,4\n10,20,30,40"
        yz = np.zeros((3, 5), dtype=float)
        result = set_friction_vegetation_values(
            yz, None, None, None, None, None, vegetation_table
        )
        expected = np.zeros_like(yz)
        expected[:, 3] = [8, 8000, 0]
        expected[:, 4] = [3, 30, 0]
        np.testing.assert_array_equal(result, expected)

    def test_with_cross_section_vegetation_values(self):
        vegetation_stem_density = 2.0
        vegetation_stem_diameter = 3.0
        vegetation_drag_coefficient = 0.5
        vegetation_height = 1.0
        yz = np.zeros((2, 5), dtype=float)
        result = set_friction_vegetation_values(
            yz,
            None,
            vegetation_stem_density,
            vegetation_stem_diameter,
            vegetation_height,
            vegetation_drag_coefficient,
            None,
        )
        expected = [[0, 0, 0, 3, 1], [0, 0, 0, 0, 0]]
        np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize(
    "width,height,match",
    [
        (
            [0, 0.5, 0],
            [0, 1, 0],
            r".*of closed profiles must have at least 4 coordinates.*",
        ),
        ([0, 0.5], [0, 1], r".*of open profiles must have at least 3 coordinates.*"),
        ([0, 0.5, 1.0], [0, 1, -1], r".*cannot have negative height coordinate.*"),
        (
            [0, 0.5, 1.0],
            [1, 2, 1],
            r".*must have at least one height coordinate at 0.*",
        ),
        (
            [0, 0.5, 1.0, 0.5],
            [0, 1, 2, 3],
            r".*should be closed or have increasing widths.*",
        ),
    ],
)
def test_tabulate_yz_err(width, height, match):
    with pytest.raises(SchematisationError, match=match):
        tabulate_yz("my-shape", width, height)


class TestCrossSectionDefinitionGetUnique:
    def test_for_closed_rectangle(self):
        csd_in = CrossSectionDefinitions(
            id=[100, 200, 300, 400],
            shape=[
                constants.CrossSectionShape.CLOSED_RECTANGLE.value,
                constants.CrossSectionShape.CLOSED_RECTANGLE.value,
                constants.CrossSectionShape.CLOSED_RECTANGLE.value,
                constants.CrossSectionShape.CLOSED_RECTANGLE.value,
            ],
            width=[1, 3, 1, 3],
            height=[1, 2, 1, 2],
            cross_section_table=["foo", "bar", "foo", "bar"],
            origin_table=["pipe", "weir", "pipe", "pipe"],
            origin_id=[10, 20, 30, 40],
        )
        unique_definition, _ = csd_in.get_unique()
        np.testing.assert_array_equal(unique_definition.id, [0, 1])
        np.testing.assert_array_equal(unique_definition.width, [1, 3])
        np.testing.assert_array_equal(unique_definition.height, [1, 2])

    @pytest.mark.parametrize(
        "shape",
        [
            constants.CrossSectionShape.RECTANGLE.value,
            constants.CrossSectionShape.CIRCLE.value,
            constants.CrossSectionShape.EGG.value,
            constants.CrossSectionShape.INVERTED_EGG.value,
        ],
    )
    def test_for_tabulated_shape(self, shape):
        csd_in = CrossSectionDefinitions(
            id=[100, 200],
            shape=[shape, shape],
            width=[1, 1],
            height=[1, 2],
            cross_section_table=["foo", "bar"],
            origin_table=["pipe", "weir"],
            origin_id=[10, 20],
        )
        unique_definition, _ = csd_in.get_unique()
        assert unique_definition.id == [
            0,
        ]
        assert unique_definition.width == [
            1,
        ]

    @pytest.mark.parametrize(
        "shape",
        [
            constants.CrossSectionShape.TABULATED_YZ.value,
            constants.CrossSectionShape.TABULATED_RECTANGLE.value,
            constants.CrossSectionShape.TABULATED_TRAPEZIUM.value,
        ],
    )
    def test_for_other_shapes(self, shape):
        csd_in = CrossSectionDefinitions(
            id=[100, 200, 300],
            shape=[shape, shape, shape],
            width=[1, 1, 100],
            height=[10, 21, 100],
            cross_section_table=["foo", "foo", "bar"],
            origin_table=["pipe", "weir", "weir"],
            origin_id=[10, 20, 30],
        )
        unique_definition, _ = csd_in.get_unique()
        np.testing.assert_array_equal(unique_definition.id, [0, 1])
        np.testing.assert_array_equal(
            sorted(unique_definition.cross_section_table), sorted(["foo", "bar"])
        )

    def test_mapping(self):
        csd_in = CrossSectionDefinitions(
            id=[100, 200, 300],
            shape=[
                constants.CrossSectionShape.CLOSED_RECTANGLE.value,
                constants.CrossSectionShape.CLOSED_RECTANGLE.value,
                constants.CrossSectionShape.RECTANGLE.value,
            ],
            width=[1, 1, 3],
            height=[1, 1, 10],
            origin_table=["pipe", "weir", "pipe"],
            origin_id=[10, 20, 20],
        )
        _, mapping = csd_in.get_unique()
        assert mapping == {"pipe": {10: 0, 20: 1}, "weir": {20: 0}}
