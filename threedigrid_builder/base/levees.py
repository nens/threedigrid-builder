from threedigrid_builder.base import array_of
from threedigrid_builder.constants import Material
from typing import Tuple

import pygeos


__all__ = ["Levees", "Breaches"]


class Breach:
    id: int
    levl: int  # the line id
    content_pk: int  # refers to v2_connected_pnt
    coordinates: Tuple[float, float]
    levee_id: int
    levmat: Material  # levee.material
    levbr: float  # levee.max_breach_depth


@array_of(Breach)
class Breaches:
    pass


class Levee:
    id: int
    the_geom: pygeos.Geometry
    crest_level: float
    max_breach_depth: float
    material: Material


@array_of(Levee)
class Levees:
    pass
