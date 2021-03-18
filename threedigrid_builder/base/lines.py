from .array import array_of
from threedigrid_builder.constants import CalculationType
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
    lik: int
    lim: int
    lin: int
    ds1d: float  # arclength
    line_geometries: pygeos.Geometry
    line_coords: Tuple[float, float, float, float]
    content_type: ContentType
    content_pk: int
    dpumax: float  # bottom_level
    flod: float  # obstacle height
    flou: float  # obstacle height
    cross1: int  # the id of the cross section location
    cross2: int  # the id of the cross section location
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
        """Set line_coords from the node coordinates"""
        start = nodes.coordinates[nodes.id_to_index(self.line[:, 0])]
        end = nodes.coordinates[nodes.id_to_index(self.line[:, 1])]
        self.line_coords[:, :2] = start
        self.line_coords[:, 2:] = end

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
