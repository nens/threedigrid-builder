import logging

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
from threedigrid_builder.constants import CalculationType, ContentType, Material

from .channels import Channels

__all__ = [
    "PotentialBreaches",
    "PotentialBreachPoints",
]

logger = logging.getLogger(__name__)


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

        # for each primary point, check the other once and keep track of the
        # secondary one
        dropped = []
        secondary_content_pk = np.full(len(self), fill_value=-9999, dtype=np.int32)
        for rng in np.split(too_close, np.where(np.diff(too_close) != 1)[0] + 1):
            if len(rng) == 0:
                continue
            if len(rng) > 1:
                dropped.extend(self.content_pk[rng[1:] + 1].tolist())
            secondary_content_pk[rng[0]] = self.content_pk[rng[0] + 1]

        if dropped:
            logger.warning(
                f"The following potential breaches will be ignored: " f"{dropped}."
            )

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
            nodes, lines, line_mask=lines.content_type == ContentType.TYPE_V2_CHANNEL
        )
        endpoints.reorder_by("node_id")

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
            idx = is_double[0] if len(is_double) > 0 else 0
            nodes.breach_ids[node_idx_] = [breach_id_1[idx], breach_id_2[idx]]

    @staticmethod
    def match_breach_ids_with_calculation_types(nodes: Nodes):
        """Make sure that the number of breach ids on a node matches its type.

        - max 1 breach for CONNECTED
        - max 2 breaches for DOUBLE_CONNECTED
        - no breaches otherwise

        Emits warnings through the logger.
        """
        has_a_breach = nodes.breach_ids[:, 0] != -9999
        has_2_breach = nodes.breach_ids[:, 1] != -9999
        is_double_connected = nodes.calculation_type == CalculationType.DOUBLE_CONNECTED
        is_single_connected = nodes.calculation_type == CalculationType.CONNECTED

        one_is_too_much_mask = has_a_breach & ~(
            is_single_connected | is_double_connected
        )
        one_is_too_much = np.where(one_is_too_much_mask)[0]
        if len(one_is_too_much) > 0:
            logger.warning(
                f"The following objects have potential breaches, but are not "
                f"(double) connected: {nodes.format_message(one_is_too_much)}."
            )
        two_is_too_much = np.where(
            has_2_breach & ~is_double_connected & ~one_is_too_much_mask
        )[0]
        if len(two_is_too_much) > 0:
            logger.warning(
                f"The following objects have two potential breaches at the "
                f"same position, but are not double connected: "
                f"{nodes.format_message(two_is_too_much)}."
            )
        ignored_breaches = np.concatenate(
            [
                nodes.breach_ids[one_is_too_much, :].ravel(),
                nodes.breach_ids[two_is_too_much, 1],
            ]
        )
        if len(ignored_breaches) > 0:
            ignored_breaches = np.unique(ignored_breaches)
            if ignored_breaches[0] == -9999:
                ignored_breaches = ignored_breaches[1:]
            logger.warning(
                f"The following potential breaches will be ignored: "
                f"{ignored_breaches.tolist()}."
            )
            nodes.breach_ids[one_is_too_much, :] = -9999
            nodes.breach_ids[two_is_too_much, 1] = -9999


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
