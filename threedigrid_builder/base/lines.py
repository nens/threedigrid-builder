from .array import array_of
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from typing import Tuple

import numpy as np
import pygeos


__all__ = ["Lines"]


class Line:
    id: int
    line_type: LineType  # kcu
    line: Tuple[int, int]
    ds1d: float  # arclength
    line_geometries: pygeos.Geometry
    content_type: ContentType
    content_pk: int
    dpumax: float  # bottom_level
    flod: float  # obstacle height
    flou: float  # obstacle height
    cross1: int  # the id of the cross section location
    cross2: int  # the id of the cross section location
    cross_weight: float


@array_of(Line)
class Lines:
    """Line between two calculation nodes (a.k.a. velocity point)."""

    def fix_line_geometries(self, nodes):
        """Construct line_geometries from node coordinates, where necessary"""
        to_fix = pygeos.is_missing(self.line_geometries)
        if not to_fix.any():
            return
        start = nodes.coordinates[nodes.id_to_index(self.line[to_fix, 0])]
        end = nodes.coordinates[nodes.id_to_index(self.line[to_fix, 1])]
        self.line_geometries[to_fix] = pygeos.linestrings(
            np.concatenate([start[:, np.newaxis], end[:, np.newaxis]], axis=1)
        )
