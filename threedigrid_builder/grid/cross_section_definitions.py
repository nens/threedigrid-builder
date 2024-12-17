import numpy as np
import shapely
from threedi_schema import constants

from threedigrid_builder.base import Array
from threedigrid_builder.constants import CrossSectionShape
from threedigrid_builder.exceptions import SchematisationError

__all__ = ["CrossSectionDefinitions", "CrossSections"]


TABULATED_YZ_DECIMALS = 3


class CrossSectionDefinition:
    id: int
    code: str
    shape: CrossSectionShape
    height: float
    width: float
    friction_values: str  # comma-separated list of floats
    vegetation_stem_density: float
    vegetation_stem_diameter: float
    vegetation_height: float
    vegetation_drag_coefficient: float
    cross_section_table: str  # csv table defining the shape of a tabulated shape
    cross_section_vegetation_table: str  # csv table with cross section vegetation table
    origin_table: str  # table definition originates from
    origin_id: int  # id in origin_table where definition originates from


class CrossSectionDefinitions(Array[CrossSectionDefinition]):
    def get_unique(self):
        """
        Returns a tuple of unique cross section definitions and their mapping to original cross section definitions.

        Returns:
            Tuple[CrossSectionDefinitions, Dict[str, Dict[int, int]]]: A tuple where the first element is
            a dictionary of unique cross section definitions and the second element is a dictionary mapping
            the original tables and rows where these definitions are used.
        """
        # Map attributes used to define a cross section for each shape
        cross_section_attributes = {
            constants.CrossSectionShape.CLOSED_RECTANGLE.value: ["width", "height"],
            constants.CrossSectionShape.RECTANGLE.value: ["width"],
            constants.CrossSectionShape.CIRCLE.value: ["width"],
            constants.CrossSectionShape.EGG.value: ["width"],
            constants.CrossSectionShape.TABULATED_RECTANGLE.value: [
                "cross_section_table"
            ],
            constants.CrossSectionShape.TABULATED_TRAPEZIUM.value: [
                "cross_section_table"
            ],
            constants.CrossSectionShape.TABULATED_YZ.value: ["cross_section_table"],
            constants.CrossSectionShape.INVERTED_EGG.value: ["width"],
        }
        definition_map = {name: {} for name in np.unique(self.origin_table)}
        new_csd_dict = {
            key: np.empty(0, dtype=data.dtype) for key, data in vars(self).items()
        }
        for shape in np.unique(self.shape):
            mask = self.shape == shape
            # collect cross section definition in array
            attr_arr = np.column_stack(
                [getattr(self, attr)[mask] for attr in cross_section_attributes[shape]]
            )
            # Find unique rows and add these to the new_csd_dict with new id's
            if len(cross_section_attributes[shape]) > 1:
                u_arr, u_idx = np.unique(attr_arr, axis=0, return_index=True)
            else:
                u_arr, u_idx = np.unique(attr_arr, return_index=True)
            new_id = np.arange(len(u_idx)) + len(new_csd_dict["id"])
            for key in new_csd_dict:
                if key == "id":
                    new_csd_dict[key] = np.concatenate([new_csd_dict[key], new_id])
                else:
                    new_csd_dict[key] = np.concatenate(
                        [new_csd_dict[key], getattr(self, key)[mask][u_idx]]
                    )
            # Map unique cross section definition to table and row of origin
            for i, row in enumerate(u_arr):
                for idx in np.where(attr_arr == row)[0]:
                    definition_map[self.origin_table[mask][idx]][
                        self.origin_id[mask][idx]
                    ] = new_id[i]
        return CrossSectionDefinitions(**new_csd_dict), definition_map

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
            count_yz=0,
        )
        if len(result) == 0:
            return result

        tables = []
        tables_yz = []

        for i, self_i in enumerate(idx):
            shape = self.shape[self_i]
            if shape in [
                CrossSectionShape.TABULATED_YZ,
                CrossSectionShape.TABULATED_RECTANGLE,
                CrossSectionShape.TABULATED_TRAPEZIUM,
            ]:
                width, height = _parse_tabulated(
                    self.cross_section_table[self_i], shape
                )
            else:
                width = self.width[self_i]
                height = self.height[self_i]
            (
                result.shape[i],
                result.width_1d[i],
                result.height_1d[i],
                table,
                yz,
            ) = tabulators[shape](shape, width, height)
            if table is not None:
                result.count[i] = len(table)
                tables.append(table)
            if yz is not None:
                result.count_yz[i] = len(yz)
                yz = set_friction_vegetation_values(
                    yz,
                    self.friction_values[self_i],
                    self.vegetation_stem_density[self_i],
                    self.vegetation_stem_diameter[self_i],
                    self.vegetation_height[self_i],
                    self.vegetation_drag_coefficient[self_i],
                    self.cross_section_vegetation_table[self_i],
                )
                tables_yz.append(yz)

        result.offset[:] = np.roll(np.cumsum(result.count), 1)
        result.offset[0] = 0
        result.offset_yz[:] = np.roll(np.cumsum(result.count_yz), 1)
        result.offset_yz[0] = 0

        if len(tables) > 0:
            result.tables = np.concatenate(tables, axis=0)
        else:
            result.tables = np.empty((0, 2))

        if len(tables_yz) > 0:
            result.tables_yz = np.concatenate(tables_yz, axis=0)
        else:
            result.tables_yz = np.empty((0, 5))

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
    offset_yz: int
    count_yz: int
    # tables: Tuple[float, float] has different length so is specified on CrossSections


class CrossSections(Array[CrossSection]):
    tables = None
    tables_yz = None


def tabulate_builtin(shape, width: float, height: float):
    """Tabulate built-in shapes (rectangle, circle)

    Built-in geometries only require a width to be fully specified.

    Args:
        shape (CrossSectionShape): returned as is
        width (float): fully specifies the shape
        height (float): ignored

    Returns:
        tuple:  shape, width_1d (float), None, None, None
    """
    return shape, width, None, None, None


def tabulate_egg(shape, width: float, height: float):
    """Tabulate the egg shape.

    Args:
        shape (CrossSectionShape): ignored
        width (float): the width of the egg, defines height and increment
        height (float): ignored; height is set to 1.5 * width

    Returns:
        tuple:  TABULATED_TRAPEZIUM, width_1d (float),
                height_1d (float), table (ndarray of shape (M, 2)), None
    """
    NUM_INCREMENTS = 16

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
    return CrossSectionShape.TABULATED_TRAPEZIUM, width, height, table, None


def tabulate_inverted_egg(shape, width, height):
    """Tabulate the egg shape, upside down.

    See tabulate_egg.
    """
    type_, width, height, table, _ = tabulate_egg(shape, width, height)
    table[:, 1] = table[::-1, 1]
    return type_, width, height, table, None


def tabulate_closed_rectangle(shape, width: float, height: float):
    """Tabulate the closed rectangle shape.

    Args:
        shape (CrossSectionShape): ignored
        width (float): the width of the rectangle
        height (float): the height of the rectangle

    Returns:
        tuple:  TABULATED_RECTANGLE, width_1d (float),
                height (float), table (ndarray of shape (M, 2)), None
    """
    table = np.array([[0.0, width], [height, 0.0]], order="F")
    return CrossSectionShape.TABULATED_RECTANGLE, width, height, table, None


def _parse_tabulated(cross_section_table, shape):
    try:
        left, right = zip(
            *[
                [float(item) for item in line.split(",")]
                for line in cross_section_table.splitlines()
            ]
        )
    except ValueError:
        raise SchematisationError(
            f"Unable to parse cross section definition table "
            f"(got: '{cross_section_table}')."
        )
    if len(left) == 0:
        raise SchematisationError(
            f"Cross section definitions of tabulated or profile type must have at least one "
            f"height element (got: {cross_section_table})."
        )
    if shape == CrossSectionShape.TABULATED_YZ:
        # for tabulated_yz, cross section table is y,z
        return np.array(left), np.array(right)
    else:
        # for other tabulated, cross seciton table is height, width
        return np.array(right), np.array(left)


def tabulate_tabulated(shape, width, height):
    """Tabulate the tabulated shapes.

    Args:
        shape (CrossSectionShape): returned as is
        width (str): space-separated widths
        height (str): space-separated heights

    Returns:
        tuple:  shape, width_1d (float),
                height_1d (float), table (ndarray of shape (M, 2)), None
    """
    if len(height) > 1 and np.any(np.diff(height) < 0.0):
        raise SchematisationError(
            f"Cross section definitions of tabulated type must have increasing heights "
            f"(got: {height})."
        )

    return shape, np.max(width), np.max(height), np.array([height, width]).T, None


def tabulate_yz(shape, width, height):
    """Tabulate an (open or closed) YZ profile

    Args:
        shape: ignored
        width (str): space-separated horizontal (Y or X) coordinates
        height (str): space-separated Z coordinates

    Returns:
        tuple:  shape, width_1d (float),
                height_1d (float), table (ndarray of shape (M, 2)), yz (ndarray of shape (M, 4))
    """
    ys = np.array(width)
    zs = np.array(height)
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
        yz = None
        shape_return = CrossSectionShape.TABULATED_TRAPEZIUM
    else:
        yz = np.zeros((len(ys), 5), dtype=float)
        yz[:, 0] = ys
        yz[:, 1] = zs
        shape_return = CrossSectionShape.TABULATED_YZ

    # Adapt non-unique height coordinates. Why?
    # Because if a segment of the profile is exactly horizontal, we need 2 widths
    seen = set()
    eps = 1 / (10 ** (TABULATED_YZ_DECIMALS + 1))
    for i, x in enumerate(zs):
        while x in seen:
            x += eps
        seen.add(x)
        zs[i] = x

    # For open profiles where the left or right sides are not the maximum Z,
    # we add a vertical line at each end such that the left and right sides are
    # now the maximum Z. This ensures that the profile is not inverted by
    # shapely when calculating the widths.
    if not is_closed:
        highest_z = np.max(zs)
        if zs[0] < highest_z:
            zs = np.concatenate(([highest_z], zs))
            ys = np.concatenate(([ys[0]], ys))
        if zs[-1] < highest_z:
            zs = np.concatenate((zs, [highest_z]))
            ys = np.concatenate((ys, [ys[-1]]))

    # shapely will automatically close an open profile
    profile = shapely.make_valid(shapely.polygons(np.array([ys, zs]).T))

    # take the length of the intersection with a horizontal line at each Z
    heights = np.unique(shapely.get_coordinates(profile)[:, 1])
    table = np.empty((len(heights), 2), dtype=float)
    y_min, y_max = ys.min(), ys.max()
    for i, height in enumerate(heights):
        line = shapely.linestrings([[y_min, height], [y_max, height]])
        cross_section_line = shapely.intersection(profile, line)
        width = shapely.length(cross_section_line)
        table[i, 0] = height
        table[i, 1] = width

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
    return shape_return, width_1d, height_1d, table, yz


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


def set_friction_vegetation_values(
    yz,
    friction_values,
    vegetation_stem_density,
    vegetation_stem_diameter,
    vegetation_height,
    vegetation_drag_coefficient,
    cross_section_vegetation_table,
):
    """Convert friction and vegetation properties from list into arrays, if available,
    and add to yz"""
    if friction_values:
        fric = np.array([float(x) for x in friction_values.split(",")])
        yz[:-1, 2] = fric
    if cross_section_vegetation_table:
        parsed = np.array(
            [
                np.fromstring(row, sep=",")
                for row in cross_section_vegetation_table.splitlines()
            ]
        )
        vegetation_stem_density = parsed[:, 0]
        vegetation_stem_diameter = parsed[:, 1]
        vegetation_height = parsed[:, 2]
        vegetation_drag_coefficient = parsed[:, 3]
    if vegetation_drag_coefficient is not None:
        yz[:-1, 3] = (
            vegetation_stem_density
            * vegetation_stem_diameter
            * vegetation_drag_coefficient
        )
        yz[:-1, 4] = vegetation_height
    return yz
