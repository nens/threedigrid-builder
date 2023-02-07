import logging

import numpy as np
import shapely

from threedigrid_builder.base import (
    Array,
    LineHalfs,
    Lines,
    Nodes,
    PointsOnLine,
    search,
)
from threedigrid_builder.constants import CalculationType, ContentType, Material
from threedigrid_builder.exceptions import SchematisationError

from .channels import Channels

__all__ = [
    "PotentialBreaches",
    "PotentialBreachPoints",
]

logger = logging.getLogger(__name__)


class PotentialBreach:
    id: int
    the_geom: shapely.Geometry
    code: str
    display_name: str
    exchange_level: float
    levee_material: Material
    maximum_breach_depth: float
    channel_id: int
    line_id: int  # only used for the output
    content_pk: int  # only used for the output


class PotentialBreachPoints(PointsOnLine):
    def merge(self) -> "PotentialBreachPoints":
        """Merge breach points that are on the same position.

        - Breaches are merged into channel start / end if they are exactly on it.
        - Breaches are merged with one another if they are exactly on the same place.
        """
        s1d = self.s1d[:]

        eps = 1e-4

        # snap breaches to channel starts/ends (to fix numerical imprecisions)
        s1d[s1d < eps] = 0.0
        lengths = self.linestrings.length[self.linestring_idx]
        mask = lengths - s1d < eps
        s1d[mask] = lengths[mask]

        # compute interdistances and a list of points that are the same
        dists = np.diff(s1d)
        same_channel = np.diff(self.linestring_idx) == 0
        same = np.where((dists < eps) & same_channel)[0]
        to_delete = same + 1  # deleting these leaves only 'primary points'

        # for each primary point keep track of the secondary one
        secondary_content_pk = np.full(len(self), fill_value=-9999, dtype=np.int32)
        for rng in np.split(same, np.where(np.diff(same) != 1)[0] + 1):
            if len(rng) == 0:
                continue
            secondary_content_pk[rng[0]] = self.content_pk[rng[0] + 1]

        # get the second one of a too close pair
        return self.__class__(
            linestrings=self.linestrings,
            id=np.delete(self.id, to_delete),
            content_pk=np.delete(self.content_pk, to_delete),
            secondary_content_pk=np.delete(secondary_content_pk, to_delete),
            linestring_idx=np.delete(self.linestring_idx, to_delete),
            s1d=np.delete(s1d, to_delete),
        )

    def find_for_line_halfs(self, line_halfs: LineHalfs):
        """Return a breach for each line_half by matching against channel id.

        It is assumed there is exactly 1 breach per line_half. This is ensured
        by PotentialBreachPoints.merge().
        """
        # per (is_start) line_half, look for a breach point (at start)
        idx = np.full(len(line_halfs), fill_value=-9999, dtype=np.int32)
        idx[line_halfs.is_start] = search(
            self.linestring_id,  # channel id
            line_halfs.content_pk[line_halfs.is_start],  # lines.content_pk
            mask=self.at_start,
            check_exists=False,
        )
        # per (is_end) line_half, look for a breach point (at end)
        idx[line_halfs.is_end] = search(
            self.linestring_id,
            line_halfs.content_pk[line_halfs.is_end],
            mask=self.at_end,
            check_exists=False,
        )
        has_breach_point = idx != -9999
        idx = idx[has_breach_point]
        return line_halfs.node_id[has_breach_point], idx

    def assign_to_connection_nodes(self, nodes: Nodes, lines: Lines):
        """Per connection node, assign max two potential breach ids.

        The priority is as follows:
        - Take the breach points of the first channel that has 2 breach
          points at the node.
        - If there are no double breach points: take the breach points of
          the first channel.
        """
        # disassemble lines into Channel - Connection Node line_halfs
        line_halfs = LineHalfs.for_connection_nodes(
            nodes, lines, line_mask=lines.content_type == ContentType.TYPE_V2_CHANNEL
        )
        line_halfs.reorder_by("node_id")

        # per line_half, match a breach point by their channel ids
        node_ids, breach_point_idx = self.find_for_line_halfs(line_halfs)
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
        """Check if the number of breach ids on a point matches its calculation type.

        - max 1 breach for CONNECTED
        - max 2 breaches for DOUBLE_CONNECTED
        - no breaches otherwise
        """
        has_a_breach = nodes.breach_ids[:, 0] != -9999
        has_2_breach = nodes.breach_ids[:, 1] != -9999
        is_double_connected = nodes.calculation_type == CalculationType.DOUBLE_CONNECTED
        is_single_connected = nodes.calculation_type == CalculationType.CONNECTED

        one_is_too_much = np.where(
            has_a_breach & ~(is_single_connected | is_double_connected)
        )[0]
        if len(one_is_too_much) > 0:
            raise SchematisationError(
                f"The following objects have potential breaches, but are not "
                f"(double) connected: {nodes.format_message(one_is_too_much)}."
            )
        two_is_too_much = np.where(has_2_breach & ~is_double_connected)[0]
        if len(two_is_too_much) > 0:
            raise SchematisationError(
                f"The following objects have two potential breaches at the "
                f"same position, but are not double connected: "
                f"{nodes.format_message(two_is_too_much)}."
            )


class PotentialBreaches(Array[PotentialBreach]):
    @property
    def side_1d(self):
        return shapely.get_point(self.the_geom, 0)

    @property
    def side_2d(self):
        return shapely.get_point(self.the_geom, -1)

    def project_on_channels(self, channels: Channels) -> PotentialBreachPoints:
        """Project the potential breaches on channels, yielding points on channels.

        This method also calls the 'merge' logic.
        """
        return PotentialBreachPoints.from_geometries(
            channels.linestrings,
            points=self.side_1d,
            linestring_idx=channels.id_to_index(self.channel_id, check_exists=True),
            content_pk=self.id,
        )
