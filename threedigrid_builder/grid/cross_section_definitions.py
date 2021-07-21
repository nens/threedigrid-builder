from threedigrid_builder.base import array_of
from threedigrid_builder.constants import CrossSectionShape
from threedigrid_builder.exceptions import SchematisationError

import math
import numpy as np


__all__ = ["CrossSectionDefinitions"]


class CrossSectionDefinition:
    id: int
    code: str
    shape: CrossSectionShape
    height: str  # space-separated list of floats
    width: str  # space-separated list of floats


class InternalCrossSectionDefinition:
    id: int
    code: str
    shape: CrossSectionShape
    content_pk: int
    width_1d: float
    offset: int
    count: int


@array_of(CrossSectionDefinition)
class CrossSectionDefinitions:
    def to_internal(self):
        """Convert to InternalCrossSectionDefinitions."""
        internal = InternalCrossSectionDefinitions(
            id=range(len(self.id)),
            content_pk=self.id,
            code=self.code,
        )
        offset = 0
        tables = []
        errored = []
        for i, shape in enumerate(self.shape):
            tabulator = tabulators[shape]
            try:
                internal.shape[i], internal.width_1d[i], table = tabulator(
                    shape, self.width[i], self.height[i]
                )
            except ValueError:
                errored.append(self.id[i])
            if table is not None:
                internal.count[i] = len(table)
                internal.offset[i] = offset
                offset += len(table)
                tables.append(table)

        if errored:
            raise SchematisationError(
                f"Unable to parse cross section definitions {errored}."
            )

        if len(tables) > 0:
            internal.tables = np.concatenate(tables, axis=0)
        else:
            internal.tables = np.empty((0, 2), order="F")

        return internal


@array_of(InternalCrossSectionDefinition)
class InternalCrossSectionDefinitions:
    def __init__(self, tables=None, **kwargs):
        self.tables = tables
        super().__init__(**kwargs)


def tabulate_builtin(shape, width, height):
    """Tabulate built-in shapes (rectangle, circle)

    Built-in geometries only require a width to be fully specified.

    Args:
        shape (CrossSectionShape): returned as is
        width (str): fully specifies the shape
        height (str): ignored

    Returns:
        tuple of TABULATED_TRAPEZIUM, width_1d (float), table (ndarray of shape (M, 2))
    """
    return shape, float(width), None


def tabulate_egg(shape, width, height):
    """Tabulate the egg shape.

    Args:
        shape (CrossSectionShape): ignored
        width (str): the width of the egg, defines height and increment
        height (str): ignored; height is set to 1.5 * width

    Returns:
        tuple of shape, width_1d (float), table (ndarray of shape (M, 2))
    """
    width = float(width)

    ## code below is copied from Inpy

    # width is the only constant
    height = width * 1.5
    increment = height / 15
    position = (2.0 / 3 * height) + (-height / 2)

    pre_heights = []
    calculating_height = -(height / 2)

    pre_heights.append(calculating_height)
    for x in range(1, 15):
        calculating_height += increment
        pre_heights.append(calculating_height)

    heights = []
    width_list = []
    for entry in pre_heights:
        calc_height = round(((entry * -1) + (height / 2.0)), 3)
        heights.append(calc_height)
        left = (-(math.pow(entry, 2)) + math.pow(height / 2.0, 2)) * math.pow(
            width / 2.0, 2
        )
        right = (
            math.pow((height / 2.0), 2) + 2.0 * position * entry + math.pow(position, 2)
        )
        eq = round(math.sqrt(left / right), 3) * 2
        width_list.append(eq)

    # start at 0.0
    heights.append(0.0)
    width_list.append(0.0)

    # flip list contents
    heights.reverse()
    width_list.reverse()

    table = np.array([heights, width_list]).T
    width_1d = np.max(table[:, 1])

    return CrossSectionShape.TABULATED_TRAPEZIUM, width_1d, table


def tabulate_closed_rectangle(shape, width, height):
    """Tabulate the closed rectangle shape.

    Args:
        shape (CrossSectionShape): ignored
        width (str): the width of the rectangle
        height (str): the height of the rectangle

    Returns:
        tuple of TABULATED_RECTANGLE, width_1d (float), table (ndarray of shape (M, 2))
    """
    width = float(width)
    height = float(height)
    table = np.array([[0.0, width], [height, 0.0]], order="F")
    return CrossSectionShape.TABULATED_RECTANGLE, width, table


def tabulate_tabulated(shape, width, height):
    """Tabulate the tabulated shapes.

    Args:
        shape (CrossSectionShape): returned as is
        width (str): space-separated widths
        height (str): space-separated heights

    Returns:
        tuple of shape, width_1d (float), table (ndarray of shape (M, 2))
    """
    table = np.array(
        [
            [float(x) for x in height.split(" ")],
            [float(x) for x in width.split(" ")],
        ]
    ).T
    width_1d = np.max(table[:, 1])
    return shape, width_1d, table


tabulators = {
    CrossSectionShape.CLOSED_RECTANGLE: tabulate_closed_rectangle,
    CrossSectionShape.RECTANGLE: tabulate_builtin,
    CrossSectionShape.CIRCLE: tabulate_builtin,
    CrossSectionShape.EGG: tabulate_egg,
    CrossSectionShape.TABULATED_RECTANGLE: tabulate_tabulated,
    CrossSectionShape.TABULATED_TRAPEZIUM: tabulate_tabulated,
}
