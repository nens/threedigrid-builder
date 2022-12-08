import numpy as np
import pygeos

from threedigrid_builder.base import Array
from threedigrid_builder.constants import CrossSectionShape
from threedigrid_builder.exceptions import SchematisationError

__all__ = ["CrossSectionDefinitions", "CrossSections"]


TABULATED_YZ_DECIMALS = 3


class CrossSectionDefinition:
    id: int
    code: str
    shape: CrossSectionShape
    height: str  # space-separated list of floats
    width: str  # space-separated list of floats


class CrossSectionDefinitions(Array[CrossSectionDefinition]):
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
        if len(result) == 0:
            return result

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


class CrossSections(Array[CrossSection]):
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
    p = (h**2 - (x**2)) * w**2
    q = h**2 + 2 * a * x + a**2
    widths = np.sqrt(p / q) * 2

    table = np.array([heights, widths]).T
    return CrossSectionShape.TABULATED_TRAPEZIUM, width, height, table


def tabulate_inverted_egg(shape, width, height):
    """Tabulate the egg shape, upside down.

    See tabulate_egg.
    """
    type_, width, height, table = tabulate_egg(shape, width, height)
    table[:, 1] = table[::-1, 1]
    return type_, width, height, table


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


def _parse_tabulated(width, height):
    try:
        heights = np.array([float(x) for x in height.split(" ")])
        widths = np.array([float(x) for x in width.split(" ")])
    except ValueError:
        raise SchematisationError(
            f"Unable to parse cross section definition width and/or height "
            f"(got: '{width}', '{height}')."
        )
    if len(heights) == 0:
        raise SchematisationError(
            f"Cross section definitions of tabulated or profile type must have at least one "
            f"height element (got: {height})."
        )
    if len(heights) != len(widths):
        raise SchematisationError(
            f"Cross section definitions of tabulated or profile type must have equal number of "
            f"height and width elements (got: {height}, {width})."
        )
    return widths, heights


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
    widths, heights = _parse_tabulated(width, height)
    if len(heights) > 1 and np.any(np.diff(heights) < 0.0):
        raise SchematisationError(
            f"Cross section definitions of tabulated type must have increasing heights "
            f"(got: {height})."
        )

    return shape, np.max(widths), np.max(heights), np.array([heights, widths]).T


def tabulate_yz(shape, width, height):
    """Tabulate an (open or closed) YZ profile

    Args:
        shape: ignored
        width (str): space-separated horizontal (Y or X) coordinates
        height (str): space-separated Z coordinates

    Returns:
        tuple:  shape, width_1d (float),
                height_1d (float), table (ndarray of shape (M, 2))
    """
    ys, zs = _parse_tabulated(width, height)
    is_closed = ys[0] == ys[-1] and zs[0] == zs[-1]
    if is_closed and len(zs) < 4:
        raise SchematisationError(
            f"Cross section definitions of closed profiles must have at least "
            f"4 coordinates (got: {len(zs)})."
        )
    if (not is_closed) and len(zs) < 3:
        raise SchematisationError(
            f"Cross section definitions of open profiles must have at least "
            f"3 coordinates (got: {len(zs)})."
        )
    if np.any(zs < 0):
        raise SchematisationError(
            f"Cross section definitions of profiles cannot have negative "
            f"height coordinate (got: {height})."
        )
    if not np.any(zs == 0):
        raise SchematisationError(
            f"Cross section definitions of profiles must have at least one "
            f"height coordinate at 0.0 (got: {height})."
        )
    if not is_closed and np.any(np.diff(ys) < 0.0):
        raise SchematisationError(
            f"Cross section definition should be closed or have increasing widths "
            f"(got: {width})."
        )

    if is_closed:
        ys = ys[:-1]
        zs = zs[:-1]

    # Adapt non-unique height coordinates. Why?
    # Because if a segment of the profile is exactly horizontal, we need 2 widths
    seen = set()
    eps = 1 / (10 ** (TABULATED_YZ_DECIMALS + 1))
    for i, x in enumerate(zs):
        while x in seen:
            x += eps
        seen.add(x)
        zs[i] = x

    # pygeos will automatically close an open profile
    profile = pygeos.make_valid(pygeos.polygons(np.array([ys, zs]).T))

    # take the length of the intersection with a horizontal line at each Z
    heights = np.unique(pygeos.get_coordinates(profile)[:, 1])
    table = np.empty((len(heights), 2), dtype=float)
    y_min, y_max = ys.min(), ys.max()
    for i, height in enumerate(heights):
        line = pygeos.linestrings([[y_min, height], [y_max, height]])
        cross_section_line = pygeos.intersection(profile, line)
        width = pygeos.length(cross_section_line)
        table[i, 0] = height
        table[i, 1] = width

    # For open profiles, if the end coordinate is closed, that means that left
    # and right sides are not equal in Z. We take the max width in that case.
    if not is_closed and table[-1, 1] == 0.0:
        table[-1, 1] = y_max - y_min

    # Eliminate duplicates and get rid of the epsilon introduced earlier
    # NB: Calccore allows a dicontinuity like [[0, 1], [1, 2], [1, 3]]
    table = np.round(table, decimals=TABULATED_YZ_DECIMALS)
    table = table[np.sort(np.unique(table, axis=0, return_index=True)[1])]

    # Drop first elements until we have 1 0.0 height at the start.
    # NB: An strong increase of width at the dry-wet transition may give issues.
    while True:
        if len(table) <= 1 or table[1, 0] > 0:
            break
        table = table[1:]

    height_1d, width_1d = table.max(axis=0).tolist()
    return CrossSectionShape.TABULATED_TRAPEZIUM, width_1d, height_1d, table


tabulators = {
    CrossSectionShape.CLOSED_RECTANGLE: tabulate_closed_rectangle,
    CrossSectionShape.RECTANGLE: tabulate_builtin,
    CrossSectionShape.CIRCLE: tabulate_builtin,
    CrossSectionShape.EGG: tabulate_egg,
    CrossSectionShape.TABULATED_RECTANGLE: tabulate_tabulated,
    CrossSectionShape.TABULATED_TRAPEZIUM: tabulate_tabulated,
    CrossSectionShape.TABULATED_YZ: tabulate_yz,
    CrossSectionShape.INVERTED_EGG: tabulate_inverted_egg,
}
