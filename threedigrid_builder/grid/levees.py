from typing import Tuple

import numpy as np
import pygeos

from threedigrid_builder.base import Array, Nodes, PointsOnLine
from threedigrid_builder.constants import Material

from .channels import Channels
from .obstacles import Obstacles

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
    node_id: int  # temporary field, for internal usage


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

        # compute interdistances and a list of points that are 'too close'
        dists = np.diff(s1d)
        same_channel = np.diff(points.linestring_idx) == 0
        too_close = np.where((dists < tolerance) & same_channel)[0]
        to_delete = too_close + 1  # deleting these leaves only 'primary points'

        # for each primary point, compute a possible secondary one
        secondary_content_pk = np.full(len(points), fill_value=-9999, dtype=np.int32)
        for rng in np.split(too_close, np.where(np.diff(too_close) != 1)[0] + 1):
            if len(rng) == 0:
                continue
            secondary_content_pk[rng[0]] = points.content_pk[rng[0] + 1]

        # adapt the s1d to be in the middle of a pair
        s1d[too_close] = (s1d[too_close] + s1d[too_close + 1]) / 2

        # get the second one of a too close pair
        return PointsOnLine(
            points.linestrings,
            id=np.delete(points.id, to_delete),
            content_pk=np.delete(points.content_pk, to_delete),
            secondary_content_pk=np.delete(secondary_content_pk, to_delete),
            linestring_idx=np.delete(points.linestring_idx, to_delete),
            s1d=np.delete(s1d, to_delete),
        )

    def set_to_connection_nodes(self, nodes: Nodes, breach_points: PointsOnLine):
        pass


class Levee:
    id: int
    the_geom: pygeos.Geometry
    crest_level: float
    max_breach_depth: float
    material: Material


class Levees(Array[Levee]):
    def merge_into_obstacles(self, obstacles: Obstacles) -> Obstacles:
        """Merge the levees into obstacles.

        This drops the 'id' column
        """
        if len(obstacles) == 0:
            first_id = 1
        else:
            first_id = obstacles.id.max() + 1
        return obstacles + Obstacles(
            id=range(first_id, first_id + len(self)),
            the_geom=self.the_geom,
            crest_level=self.crest_level,
        )
