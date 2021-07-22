from threedigrid_builder.base import array_of
from threedigrid_builder.constants import CrossSectionShape
from threedigrid_builder.exceptions import SchematisationError

import numpy as np


__all__ = ["CrossSectionDefinitions", "InternalCrossSectionDefinitions"]


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
    tables = None


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
    # width is the only constant; height derives from width
    width = float(width)
    height = width * 1.5

    heights = np.linspace(height / 15, height, num=15, endpoint=True)

    position = height / 6
    pre_heights = (height / 2) - heights
    left = ((height / 2) ** 2 - (pre_heights ** 2)) * (width / 2) ** 2
    right = (height / 2) ** 2 + 2 * position * pre_heights + position ** 2
    widths = np.sqrt(left / right) * 2

    table = np.empty((16, 2))
    table[0] = 0.0
    table[1:, 0] = np.round(heights, 3)
    table[1:, 1] = np.round(widths, 3)
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
