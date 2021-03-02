from threedigrid_builder.base import array_of
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import CrossSectionShape
from threedigrid_builder.constants import FrictionType

import numpy as np
import pygeos


__all__ = ["CrossSectionLocations", "CrossSectionDefinitions"]


class CrossSectionLocation:
    id: int
    # code: str  unused?
    the_geom: pygeos.Geometry
    definition_id: id  # refers to CrossSectionDefinition
    channel_id: id  # refers to Channel
    reference_level: float
    bank_level: float
    friction_type: FrictionType
    friction_value: float


@array_of(CrossSectionLocation)
class CrossSectionLocations:
    def apply_to_channels(self, channels, lines):
        """Compute cross section weights for channels.

        As a side-effect, this function sorts self inplace on (channel_id,
        position on channel).

        Args:
            channels (Channels): will be changed in place
            lines (Lines): will be changed in place
        """
        if not np.in1d(channels.id, self.channel_id).all():
            missing = set(channels.id) - set(self.channel_id)
            raise ValueError(f"Channels {missing} have no cross section location set.")

        # get the indices that sort channel lines by their content_pk
        line_idx = np.lexsort((lines.internal_seq_id, lines.content_pk))
        line_idx = line_idx[lines.content_type == ContentType.TYPE_V2_CHANNEL]

        # reorder self by channel_id
        self.reorder_by(np.argsort(self.channel_id))
        _, n_locations = np.unique(self.channel_id, return_counts=True)

        # compute where a location is on each channel, cumulative over
        # all channels
        location_ds = pygeos.line_locate_point(
            np.repeat(channels.the_geom, n_locations), self.the_geom
        )
        channel_ds_cumul = np.cumsum(pygeos.length(channels.the_geom))
        channel_ds_cumul = np.roll(channel_ds_cumul, 1)
        channel_ds_cumul[0] = 0
        location_ds += np.repeat(channel_ds_cumul, n_locations)

        # reorder self again so that locations are ordered on each channel
        _location_ds_sorter = np.argsort(location_ds)
        location_ds = location_ds[_location_ds_sorter]
        self.reorder_by(_location_ds_sorter)

        # compute where a line-midpoint (velocity point) is, cumulative over
        # all channels
        cumul = np.cumsum(lines.ds1d[line_idx])
        midpoints = (cumul + np.roll(cumul, 1)) / 2
        midpoints[0] = cumul[0] / 2

        # for each midpoint, find what location comes after
        location_idx_after = np.searchsorted(location_ds, midpoints)
        location_idx_before = location_idx_after - 1

        # clip so that we do not get out-of-bounds errors
        location_idx_after[location_idx_after >= len(self)] = len(self) - 1
        location_idx_before[location_idx_before < 0] = 0

        lines.cross1[line_idx] = self.id[location_idx_before]
        lines.cross2[line_idx] = self.id[location_idx_after]

        # make cross1 and cross2 equal if they are not from the correct channel
        to_equalize = self.channel_id[location_idx_before] != lines.content_pk[line_idx]
        lines.cross1[to_equalize] = lines.cross2[to_equalize]
        to_equalize = self.channel_id[location_idx_after] != lines.content_pk[line_idx]
        lines.cross2[to_equalize] = lines.cross2[to_equalize]

        # for each cross section location, test how far it is from the channel
        # start


class CrossSectionDefinition:
    id: int
    # code: str  unused?
    shape: CrossSectionShape
    height: float
    width: float


@array_of(CrossSectionDefinition)
class CrossSectionDefinitions:
    pass
