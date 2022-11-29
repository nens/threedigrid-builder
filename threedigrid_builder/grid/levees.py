from typing import Tuple

import pygeos

from threedigrid_builder.base import Array
from threedigrid_builder.constants import Material

__all__ = ["Levees", "Breaches", "PotentialBreaches"]


class Breach:
    # Deprecated
    id: int
    levl: int  # the line id
    content_pk: int  # refers to v2_connected_pnt
    coordinates: Tuple[float, float]
    levee_id: int
    levmat: Material  # levee.material
    levbr: float  # levee.max_breach_depth


class Breaches(Array[Breach]):
    pass


class PotentialBreach:
    id: int
    the_geom: pygeos.Geometry
    code: str
    display_name: str
    levee_material: Material
    maximum_breach_depth: float
    channel_id: int
    line_id: int


class PotentialBreaches(Array[PotentialBreach]):
    @property
    def side_1d(self):
        return pygeos.get_point(self.the_geom, 0)

    @property
    def side_2d(self):
        return pygeos.get_point(self.the_geom, -1)


class Levee:
    id: int
    the_geom: pygeos.Geometry
    crest_level: float
    max_breach_depth: float
    material: Material


class Levees(Array[Levee]):
    pass
