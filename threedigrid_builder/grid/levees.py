from typing import Tuple

import numpy as np
import pygeos

from threedigrid_builder.base import Array, PointsOnLine
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
    # levee_material: Material
    # maximum_breach_depth: float
    channel_id: int
    line_id: int


class PotentialBreaches(Array[PotentialBreach]):
    @property
    def side_1d(self):
        return pygeos.get_point(self.the_geom, 0)

    @property
    def side_2d(self):
        return pygeos.get_point(self.the_geom, -1)

    def project_on_channels(
        self, channels: Channels, merge_tolerance: float
    ) -> PointsOnLine:
        """Project the potential breaches on channels, yielding points on channels.

        This method also calls the 'merge' logic.
        """
        points = PointsOnLine.from_geometries(
            channels.linestrings,
            points=self.side_1d,
            linestring_idx=channels.id_to_index(self.channel_id, check_exists=True),
            content_pk=self.id,
        )
        return self.merge(points, merge_tolerance)

    @staticmethod
    def merge(points: PointsOnLine, tolerance: float) -> PointsOnLine:
        """Merge breach points with a certain tolerance.

        - Breaches are merged into channel start / end if they are closer than
          tolerance to it.
        - Breaches are merged with one another if they are closer than tolerance
          to one another.
        - If there are more than 2 breaches too close to one another, the third
          and so forth will be dropped.
        """
        s1d = points.s1d[:]

        # snap breaches to channel starts/ends
        s1d[s1d < tolerance] = 0.0
        lengths = points.linestrings.length[points.linestring_idx]
        mask = lengths - s1d < tolerance
        s1d[mask] = lengths[mask]

        # compute interdistances
        dists = np.diff(s1d)
        same_channel = np.diff(points.linestring_idx) == 0
        too_close = np.where((dists < tolerance) & same_channel)[0]

        # disable the second one of a 'too close' pair
        mask = np.full(len(points), fill_value=True, dtype=bool)
        mask[too_close + 1] = False

        # adapt the s1d to be in the middle of a pair
        s1d[too_close] = (s1d[too_close] + s1d[too_close + 1]) / 2
        return PointsOnLine(
            points.linestrings,
            id=points.id[mask],
            content_pk=points.content_pk[mask],
            linestring_idx=points.linestring_idx[mask],
            s1d=s1d[mask],
        )


class Levee:
    id: int
    the_geom: pygeos.Geometry
    crest_level: float
    max_breach_depth: float
    material: Material


class Levees(Array[Levee]):
    pass
