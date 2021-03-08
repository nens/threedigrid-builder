from .array import array_of
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import LineType
from typing import Tuple

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
