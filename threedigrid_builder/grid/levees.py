from typing import Tuple

import numpy as np
import pygeos

from threedigrid_builder.base import (
    Array,
    Endpoints,
    Lines,
    Nodes,
    PointsOnLine,
    search,
)
from threedigrid_builder.constants import ContentType, Material

from .channels import Channels
from .obstacles import Obstacles

__all__ = [
    "Levees",
    "Breaches",
    "PotentialBreaches",
    "PotentialBreachPoints",
]


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
    exchange_level: float
    levee_material: Material
    maximum_breach_depth: float
    channel_id: int


class PotentialBreachPoints(PointsOnLine):
    def merge(self, tolerance: float) -> "PotentialBreachPoints":
        """Merge breach points with a certain tolerance.

        - Breaches are merged into channel start / end if they are closer than
          tolerance to it.
        - Breaches are merged with one another if they are closer than tolerance
          to one another.
        - If there are more than 2 breaches too close to one another, the third
          and so forth will be dropped.
        """
        s1d = self.s1d[:]

        # snap breaches to channel starts/ends
        s1d[s1d < tolerance] = 0.0
        lengths = self.linestrings.length[self.linestring_idx]
        mask = lengths - s1d < tolerance
        s1d[mask] = lengths[mask]

        # compute interdistances and a list of points that are 'too close'
        dists = np.diff(s1d)
        same_channel = np.diff(self.linestring_idx) == 0
        too_close = np.where((dists < tolerance) & same_channel)[0]
        to_delete = too_close + 1  # deleting these leaves only 'primary points'

        # for each primary point, compute a possible secondary one
        secondary_content_pk = np.full(len(self), fill_value=-9999, dtype=np.int32)
        for rng in np.split(too_close, np.where(np.diff(too_close) != 1)[0] + 1):
            if len(rng) == 0:
                continue
            secondary_content_pk[rng[0]] = self.content_pk[rng[0] + 1]

        # adapt the s1d to be in the middle of a pair
        s1d[too_close] = (s1d[too_close] + s1d[too_close + 1]) / 2

        # get the second one of a too close pair
        return self.__class__(
            linestrings=self.linestrings,
            id=np.delete(self.id, to_delete),
            content_pk=np.delete(self.content_pk, to_delete),
            secondary_content_pk=np.delete(secondary_content_pk, to_delete),
            linestring_idx=np.delete(self.linestring_idx, to_delete),
            s1d=np.delete(s1d, to_delete),
        )

    def find_for_endpoints(self, endpoints: Endpoints):
        """Return a breach for each endpoint by matching against channel id.

        It is assumed there is exactly 1 breach per endpoint. This is ensured
        by PotentialBreachPoints.merge().
        """
        # per (is_start) endpoint, look for a breach point (at start)
        idx = np.full(len(endpoints), fill_value=-9999, dtype=np.int32)
        idx[endpoints.is_start] = search(
            self.linestring_id,  # channel id
            endpoints.content_pk[endpoints.is_start],  # lines.content_pk
            mask=self.at_start,
            check_exists=False,
        )
        # per (is_end) endpoint, look for a breach point (at end)
        idx[endpoints.is_end] = search(
            self.linestring_id,
            endpoints.content_pk[endpoints.is_end],
            mask=self.at_end,
            check_exists=False,
        )
        has_breach_point = idx != -9999
        idx = idx[has_breach_point]
        return endpoints.node_id[has_breach_point], idx

    def assign_to_connection_nodes(self, nodes: Nodes, lines: Lines):
        """Per connection node, assign max two potential breach ids.

        The priority is as follows:
        - Take the breach points of the first channel that has 2 breach
          points at the node.
        - If there are no double breach points: take the breach points of
          the first channel.
        """
        # disassemble lines into Channel - Connection Node endpoints
        endpoints = Endpoints.for_connection_nodes(
            nodes, lines, line_types=[ContentType.TYPE_V2_CHANNEL]
        )

        # per endpoint, match a breach point by their channel ids
        node_ids, breach_point_idx = self.find_for_endpoints(endpoints)
        node_idx = nodes.id_to_index(node_ids)

        # iterate over nodes with a breach point and assign (if present)
        # the first double breach point
        for node_idx_ in np.unique(node_idx):
            options = breach_point_idx[node_idx == node_idx_]
            # sort by linestring (channel) id for reproducibility:
            options = options[np.argsort(self.linestring_id[options])]
            breach_id_1 = self.content_pk[options]
            breach_id_2 = self.secondary_content_pk[options]
            is_double = np.where(breach_id_2 != -9999)[0]
            if len(is_double) > 0:
                nodes.breach_ids[node_idx_] = [
                    breach_id_1[is_double[0]],
                    breach_id_2[is_double[0]],
                ]
            else:
                nodes.breach_ids[node_idx_, 0] = breach_id_1[0]


class PotentialBreaches(Array[PotentialBreach]):
    @property
    def side_1d(self):
        return pygeos.get_point(self.the_geom, 0)

    @property
    def side_2d(self):
        return pygeos.get_point(self.the_geom, -1)

    def project_on_channels(
        self, channels: Channels, merge_tolerance: float
    ) -> PotentialBreachPoints:
        """Project the potential breaches on channels, yielding points on channels.

        This method also calls the 'merge' logic.
        """
        return PotentialBreachPoints.from_geometries(
            channels.linestrings,
            points=self.side_1d,
            linestring_idx=channels.id_to_index(self.channel_id, check_exists=True),
            content_pk=self.id,
        ).merge(merge_tolerance)


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
