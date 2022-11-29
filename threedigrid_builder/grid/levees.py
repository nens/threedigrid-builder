from typing import Tuple

import pygeos

from threedigrid_builder.base import Array
from threedigrid_builder.constants import Material

from .channels import Channels

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

    def project_on_channel(self, channels: Channels):
        """Return the dist_to_start of breaches on its corresponding channel."""
        channel_idx = channels.id_to_index(self.channel_id, check_exists=True)
        dist_to_start = pygeos.line_locate_point(
            channels.the_geom[channel_idx], self.side_1d
        )
        return channel_idx, dist_to_start


class Levee:
    id: int
    the_geom: pygeos.Geometry
    crest_level: float
    max_breach_depth: float
    material: Material


class Levees(Array[Levee]):
    pass
