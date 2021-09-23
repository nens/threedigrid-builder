from .array import array_of
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from typing import Tuple

import numpy as np
import pygeos


__all__ = ["Lines"]


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
    cross1: int  # the id of the cross section definition
    cross2: int  # the id of the cross section definition
    frict_type1: int  # Friction type
    frict_type2: int  # Friction type
    frict_value1: float  # Friction type
    frict_value2: float  # Friction type
    cross_weight: float
    discharge_coefficient_positive: float
    discharge_coefficient_negative: float


@array_of(Line)
class Lines:
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

    def set_bottom_levels(self, nodes, allow_nan=False):
        """Set bottom levels (dpumax) for lines.

        The bottom level for a line is the largest of the bottom level (dmax) of the
        two nodes the line is attached to, except for weirs and orifices, where it
        equals the crest_level (should be set already to dpumax).

        Args:
            nodes (Nodes): nodes from which to take the dmax
            allow_nan (bool): whether to error if any dmax is NaN. If True, NaN is
                propagated to dpumax if any of the two dmax is NaN. Default False.
        """
        mask = ~np.isfinite(self.dpumax)  # filter out already set dpumax (weir/orifice)
        if not mask.any():
            return
        node_dmax = nodes.dmax.take(nodes.id_to_index(self.line[mask]))
        if not allow_nan and np.any(~np.isfinite(node_dmax)):
            raise ValueError("Found NaN after setting bottom levels")
        with np.errstate(invalid="ignore"):  # suppress warnings on max(nan, nan)
            self.dpumax[mask] = np.max(node_dmax, axis=1)

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
        self.ds1d[to_fix] = pygeos.length(self.line_geometries[to_fix])

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
