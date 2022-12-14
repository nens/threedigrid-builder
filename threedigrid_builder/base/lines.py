from typing import Tuple

import numpy as np
import pygeos

from threedigrid_builder.constants import CalculationType, ContentType, LineType

from .array import Array

__all__ = ["Lines", "Endpoints"]


class Line:
    id: int
    code: str
    display_name: str
    kcu: LineType  # calculation type of the line
    line: Tuple[int, int]
    cross_pix_coords: Tuple[int, int, int, int]
    lik: int  # quadtree grid coordinate z of the line start
    lim: int  # quadtree grid coordinate x of the line start
    lin: int  # quadtree grid coordinate y of the line start
    s1d: float  # position (arclength) along a 1D element
    ds1d_half: float  # position (arclength) along a single line (between two calc nodes.)
    ds1d: float  # arclength
    line_geometries: pygeos.Geometry
    line_coords: Tuple[float, float, float, float]
    content_type: ContentType
    content_pk: int
    dpumax: float  # bottom_level at the velocity point
    invert_level_start_point: float  # bottom level at line start
    invert_level_end_point: float  # bottom level at line end
    flod: float  # obstacle height
    flou: float  # obstacle height
    cross_id1: int  # the id of the cross section definition
    cross_id2: int  # the id of the cross section definition
    frict_type1: int  # Friction type
    frict_type2: int  # Friction type
    frict_value1: float  # Friction type
    frict_value2: float  # Friction type
    cross_weight: float
    discharge_coefficient_positive: float
    discharge_coefficient_negative: float
    is_1d_boundary: int  # internal flag
    zoom_category: int
    connection_node_start_id: int
    connection_node_end_id: int
    crest_level: float
    crest_type: CalculationType
    dist_calc_points: float
    material: int
    sewerage: int
    sewerage_type: int
    windshieldings: Tuple[float, float, float, float, float, float, float, float]


class Lines(Array[Line]):
    """Line between two calculation nodes (a.k.a. velocity point)."""

    def set_discharge_coefficients(self):
        """Set discharge coefficients to 1.0 where unset."""
        for arr in (
            self.discharge_coefficient_positive,
            self.discharge_coefficient_negative,
        ):
            arr[np.isnan(arr)] = 1.0

    def set_line_coords(self, nodes):
        """Set line_coords from the node coordinates where necessary"""
        to_fix = np.isnan(self.line_coords).any(axis=1)
        start = nodes.coordinates[nodes.id_to_index(self.line[to_fix, 0])]
        end = nodes.coordinates[nodes.id_to_index(self.line[to_fix, 1])]
        self.line_coords[to_fix, :2] = start
        self.line_coords[to_fix, 2:] = end

    def fix_line_geometries(self):
        """Construct line_geometries from line_coords, where necessary"""
        to_fix = pygeos.is_missing(self.line_geometries)
        if not to_fix.any():
            return
        if np.any(~np.isfinite(self.line_coords[to_fix])):
            raise ValueError("No line coords available")
        self.line_geometries[to_fix] = pygeos.linestrings(
            self.line_coords[to_fix].reshape(-1, 2, 2)
        )

    def fix_ds1d(self):
        """Construct line distance (ds1d) from line_geometries, where necessary"""
        to_fix = np.isnan(self.ds1d)
        if not to_fix.any():
            return
        if np.any(pygeos.is_missing(self.line_geometries[to_fix])):
            raise ValueError("No line geometries available")
        self.ds1d[to_fix] = pygeos.length(self.line_geometries[to_fix])
        self.ds1d_half[to_fix] = 0.5 * self.ds1d[to_fix]

    def sort_by_nodes(self, node_ids):
        """Order selected lines by node id, ascending.

        Only lines with given node ids are ordered and then on then id that is supplied.
        This operation only makes sense if every node id occurs once in lines.line.
        """
        line_idx, start_or_end = np.where(np.isin(self.line, node_ids))
        new_line_idx = np.arange(len(self))
        sorter = np.argsort(self.line[line_idx, start_or_end])
        new_line_idx[line_idx] = new_line_idx[line_idx][sorter]
        self.reorder(new_line_idx)

    def set_2d_crest_levels(self, crest_levels, where):
        """Set obstacle/levee crest levels to 2D lines

        Sets: flod, flou, kcu inplace (for 2D lines)
        """
        has_crest_level = where[np.isfinite(crest_levels)]
        is_2d_u = self.kcu[has_crest_level] == LineType.LINE_2D_U
        is_2d_v = self.kcu[has_crest_level] == LineType.LINE_2D_V
        self.kcu[has_crest_level[is_2d_u]] = LineType.LINE_2D_OBSTACLE_U
        self.kcu[has_crest_level[is_2d_v]] = LineType.LINE_2D_OBSTACLE_V
        self.flod[where] = crest_levels
        self.flou[where] = crest_levels

    def as_endpoints(self, where=None) -> "Endpoints":
        if where is None:
            where = np.arange(len(self), dtype=int)
        elif where.dtype == bool:
            where = np.where(where)[0]

        n = len(where)

        result = Endpoints(
            lines=self,
            id=range(n * 2),
            line_idx=np.repeat(where, 2),
            node_id=self.line[where].ravel(),
            is_start=np.tile([True, False], n),
        )
        result.reorder_by("node_id")
        return result


class Endpoint:
    id: int
    line_idx: int
    is_start: bool
    node_id: int


class Endpoints(Array[Endpoint]):
    """Each line has 2 endpoints."""

    scalars = ("lines",)

    @property
    def is_end(self):
        return ~self.is_start

    @property
    def invert_level(self):
        return np.where(
            self.is_start, self.invert_level_start_point, self.invert_level_end_point
        )

    def __getattr__(self, name):
        return getattr(self.lines, name)[self.line_idx]

    def filter_by_node_id(self, node_ids) -> "Endpoints":
        return self[np.isin(self.node_id, node_ids)]

    def reduce_per_node(self, ufunc, values):
        """Compute the per-node ufunc reduction of 'values'

        For computing the maximum, use 'np.fmax'. For computing the sum, use 'np.add'.

        Returns node ids and corresponding maximum values.
        """
        if len(values) != len(self):
            raise ValueError("values must have the same length as self")
        if len(values) == 0:
            return np.empty(0, dtype=int), np.empty(0, dtype=float)
        diff = np.diff(self.node_id)
        if np.any(diff < 0):
            raise ValueError("node_id must be sorted")
        indices = np.concatenate([[0], np.where(diff > 0)[0] + 1])
        result = ufunc.reduceat(values, indices)
        return self.node_id[indices], result
