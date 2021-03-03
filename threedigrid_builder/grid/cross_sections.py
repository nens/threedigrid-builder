from threedigrid_builder.base import array_of
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import CrossSectionShape
from threedigrid_builder.constants import FrictionType

import numpy as np
import pygeos


__all__ = ["CrossSectionLocations", "CrossSectionDefinitions"]


class CrossSectionLocation:
    id: int
    code: str
    the_geom: pygeos.Geometry
    definition_id: id  # refers to CrossSectionDefinition
    channel_id: id  # refers to Channel
    reference_level: float
    bank_level: float
    friction_type: FrictionType
    friction_value: float


@array_of(CrossSectionLocation)
class CrossSectionLocations:
    pass


class CrossSectionDefinition:
    id: int
    code: str
    shape: CrossSectionShape
    height: float
    width: float


@array_of(CrossSectionDefinition)
class CrossSectionDefinitions:
    pass


def compute_weights(lines, cs, channels):
    """Compute cross section weights for channel lines.

    Args:
        lines (Lines): Lines for which to compute cross1, cross2, and
            cross_weight. The attributes will be changed in place. Only
            lines of type Channel will be included. Those lines MUST be
            ordered by channel_id and then by position on channel. If this
            is not the case, the computed weights will be bogus.
            Required attributes: content_type, content_pk, and ds1d.
        cs (CrossSectionLocations)
        channels (Channels): Used to lookup the channel geometry

    See also:
        Channels.get_lines: computes the lines in the correct order
    """
    if not np.in1d(channels.id, cs.channel_id).all():
        missing = set(channels.id) - set(cs.channel_id)
        raise ValueError(f"Channels {missing} have no cross section location set.")

    # Mask the lines to only the Channel lines
    line_mask = lines.content_type == ContentType.TYPE_V2_CHANNEL

    # Make an array that sorts the cross section (cs) locations by channel_id
    cs_sorter = np.argsort(cs.channel_id)

    # Count how many cs locations we have per channel
    _, cs_per_channel = np.unique(cs.channel_id, return_counts=True)

    # Locate the position measured along the channel (ds) of the CS locations
    cs_ds = pygeos.line_locate_point(
        np.repeat(channels.the_geom, cs_per_channel), cs.the_geom[cs_sorter]
    )

    # To each cs_ds, add the length of all channels before it
    _cumulative = np.cumsum(pygeos.length(channels.the_geom))
    _cumulative = np.roll(_cumulative, 1)
    _cumulative[0] = 0
    cs_ds += np.repeat(_cumulative, cs_per_channel)

    # Make cs_ds monotonically increasing and update cs_sorter to match it
    _cs_ds_sorter = np.argsort(cs_ds)
    cs_ds = cs_ds[_cs_ds_sorter]
    cs_sorter = cs_sorter[_cs_ds_sorter]

    # Compute the ds of the line midpoints, cumulative over all channels
    _cumulative = np.cumsum(lines.ds1d[line_mask])
    midpoint_ds = (_cumulative + np.roll(_cumulative, 1)) / 2
    midpoint_ds[0] = _cumulative[0] / 2

    # Find what CS location comes after and before each midpoint
    cs_idx_2 = np.searchsorted(cs_ds, midpoint_ds)
    cs_idx_1 = cs_idx_2 - 1
    out_of_bounds_1 = cs_idx_1 < 0
    out_of_bounds_2 = cs_idx_2 >= len(cs)
    cs_idx_1[out_of_bounds_1] = 0
    cs_idx_2[out_of_bounds_2] = len(cs) - 1

    # Extrapolation: if a midpoint does not have a cross section location
    # to its left (cs_idx_1 is invalid), assign to it the 2 cross_idx
    # to its right. And vice versa for cs_idx_2 that is invalid. Note
    # that because every channel has at least 1 crosssection, one the two
    # must always be valid. This is not checked here.
    # Equalization: if a newly assigned crosssection location would not
    # belong to the correct channel, we have only 1 crossection location
    # in the channel so we assign the same cross_idx to both 1 and 2.

    # Create two matching channel id arrays
    cs_ch_id = cs.channel_id[cs_sorter]
    line_ch_id = lines.content_pk[line_mask]

    extrap_left = (cs_ch_id[cs_idx_1] != line_ch_id) | out_of_bounds_1
    cs_idx_1[extrap_left] = cs_idx_2[extrap_left]
    # additional filtering to prevent cs_idx_2 being out of bounds
    extrap_left[extrap_left] &= cs_idx_2[extrap_left] < len(cs) - 1
    extrap_left[extrap_left] &= (
        cs_ch_id[cs_idx_2[extrap_left] + 1] == line_ch_id[extrap_left]
    )
    cs_idx_2[extrap_left] += 1

    extrap_right = (cs_ch_id[cs_idx_2] != line_ch_id) | out_of_bounds_2
    cs_idx_2[extrap_right] = cs_idx_1[extrap_right]
    # additional filtering to prevent cs_idx_1 being out of bounds
    extrap_right[extrap_right] &= cs_idx_1[extrap_right] > 0
    extrap_right[extrap_right] &= (
        cs_ch_id[cs_idx_1[extrap_right] - 1] == line_ch_id[extrap_right]
    )
    cs_idx_1[extrap_right] -= 1

    # map index to id and set on the (filtered) lines
    lines.cross1[line_mask] = cs.id[cs_sorter][cs_idx_1]
    lines.cross2[line_mask] = cs.id[cs_sorter][cs_idx_2]

    # compute the weights. for each line, we have 3 times ds:
    # 1. the ds of the CrossSectionLocation before it
    # 2. the ds of the CrossSectionLocation after it
    # 3. the ds of the midpoint (velocity point)
    ds_1 = cs_ds[cs_idx_1]
    ds_2 = cs_ds[cs_idx_2]
    with np.errstate(divide="ignore", invalid="ignore"):
        # this transforms 1 / 0 to inf and 0 / 0 to nan without warning
        weights = (ds_2 - midpoint_ds) / (ds_2 - ds_1)
    weights[~np.isfinite(weights)] = 1.0

    lines.cross_weight[line_mask] = weights
