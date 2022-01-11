from threedigrid_builder.base import array_of
from threedigrid_builder.constants import CrossSectionShape
from threedigrid_builder.exceptions import SchematisationError

import numpy as np


__all__ = ["CrossSectionDefinitions", "CrossSections"]


class CrossSectionDefinition:
    id: int
    code: str
    shape: CrossSectionShape
    height: str  # space-separated list of floats
    width: str  # space-separated list of floats


@array_of(CrossSectionDefinition)
class CrossSectionDefinitions:
    def convert(self, ids):
        """Convert to CrossSections.

        Args:
            ids (ndarray of int): A list of cross section definition ids to convert.
        """
        try:
            idx = self.id_to_index(ids, check_exists=True)
        except KeyError as e:
            if len(e.values) > 10:
                id_msg = str(e.values[:10].tolist()).replace("]", ", ...]")
            else:
                id_msg = str(e.values.tolist())
            raise SchematisationError(
                f"One or more objects refer to non-existing cross section definitions "
                f"{id_msg}."
            )

        result = CrossSections(
            id=range(len(idx)),
            content_pk=ids,
            code=self.code[idx],
            count=0,
        )
        offset = 0
        tables = []

        # Numpy array views on width/height based on idx
        width_idx = self.width[idx]
        height_idx = self.height[idx]

        for i, shape in enumerate(self.shape[idx]):
            tabulator = tabulators[shape]
            result.shape[i], result.width_1d[i], result.height_1d[i], table = tabulator(
                shape, width_idx[i], height_idx[i]
            )
            if table is not None:
                result.count[i] = len(table)
                result.offset[i] = offset
                offset += len(table)
                tables.append(table)

        result.offset[:] = np.roll(np.cumsum(result.count), 1)
        result.offset[0] = 0

        if len(tables) > 0:
            result.tables = np.concatenate(tables, axis=0)
        else:
            result.tables = np.empty((0, 2))

        return result


class CrossSection:
    id: int
    code: str
    shape: CrossSectionShape
    content_pk: int
    width_1d: float
    height_1d: float
    offset: int
    count: int
    # tables: Tuple[float, float] has different length so is specified on CrossSections


@array_of(CrossSection)
class CrossSections:
    tables = None


def tabulate_builtin(shape, width, height):
    """Tabulate built-in shapes (rectangle, circle)

    Built-in geometries only require a width to be fully specified.

    Args:
        shape (CrossSectionShape): returned as is
        width (str): fully specifies the shape
        height (str): ignored

    Returns:
        tuple:  shape, width_1d (float), None, None
    """
    try:
        width = float(width)
    except ValueError:
        raise SchematisationError(
            f"Unable to parse cross section definition width (got: '{width}')."
        )

    return shape, width, None, None


def tabulate_egg(shape, width, height):
    """Tabulate the egg shape.

    Args:
        shape (CrossSectionShape): ignored
        width (str): the width of the egg, defines height and increment
        height (str): ignored; height is set to 1.5 * width

    Returns:
        tuple:  TABULATED_TRAPEZIUM, width_1d (float),
                height_1d (float), table (ndarray of shape (M, 2))
    """
    NUM_INCREMENTS = 16

    # width is the only constant; height derives from width
    try:
        width = float(width)
    except ValueError:
        raise SchematisationError(
            f"Unable to parse cross section definition width (got: '{width}')."
        )

    height = width * 1.5
    position = height / 3.0  # some parameter for the 'egg' curve
    heights = np.linspace(0, height, num=NUM_INCREMENTS, endpoint=True)

    # here comes the formulation of the "egg" curve, simplified from inpy
    h = height / 2
    w = width / 2
    a = position / 2
    x = h - heights
    p = (h ** 2 - (x ** 2)) * w ** 2
    q = h ** 2 + 2 * a * x + a ** 2
    widths = np.sqrt(p / q) * 2

    table = np.array([heights, widths]).T
    return CrossSectionShape.TABULATED_TRAPEZIUM, width, height, table


def tabulate_closed_rectangle(shape, width, height):
    """Tabulate the closed rectangle shape.

    Args:
        shape (CrossSectionShape): ignored
        width (str): the width of the rectangle
        height (str): the height of the rectangle

    Returns:
        tuple:  TABULATED_RECTANGLE, width_1d (float),
                height (float), table (ndarray of shape (M, 2))
    """
    try:
        width = float(width)
        height = float(height)
    except ValueError:
        raise SchematisationError(
            f"Unable to parse cross section definition width and/or height "
            f"(got: '{width}', '{height}')."
        )
    table = np.array([[0.0, width], [height, 0.0]], order="F")
    return CrossSectionShape.TABULATED_RECTANGLE, width, height, table


def tabulate_tabulated(shape, width, height):
    """Tabulate the tabulated shapes.

    Args:
        shape (CrossSectionShape): returned as is
        width (str): space-separated widths
        height (str): space-separated heights

    Returns:
        tuple:  shape, width_1d (float),
                height_1d (float), table (ndarray of shape (M, 2))
    """
    try:
        heights = [float(x) for x in height.split(" ")]
        widths = [float(x) for x in width.split(" ")]
    except ValueError:
        raise SchematisationError(
            f"Unable to parse cross section definition width and/or height "
            f"(got: '{width}', '{height}')."
        )
    if len(heights) == 0:
        raise SchematisationError(
            f"Cross section definitions of tabulated type must have at least one "
            f"height element (got: {height})."
        )
    if len(heights) != len(widths):
        raise SchematisationError(
            f"Cross section definitions of tabulated type must have equal number of "
            f"height and width elements (got: {height}, {width})."
        )
    if len(heights) > 1 and np.any(np.diff(heights) < 0.0):
        raise SchematisationError(
            f"Cross section definitions of tabulated type must have increasing heights "
            f"(got: {height})."
        )

    return shape, np.max(widths), np.max(heights), np.array([heights, widths]).T


tabulators = {
    CrossSectionShape.CLOSED_RECTANGLE: tabulate_closed_rectangle,
    CrossSectionShape.RECTANGLE: tabulate_builtin,
    CrossSectionShape.CIRCLE: tabulate_builtin,
    CrossSectionShape.EGG: tabulate_egg,
    CrossSectionShape.TABULATED_RECTANGLE: tabulate_tabulated,
    CrossSectionShape.TABULATED_TRAPEZIUM: tabulate_tabulated,
}
