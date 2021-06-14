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


def compute_weights(channel_id, ds, cs, channels):
    """Compute cross section weights for points on channels.

    Points on channels are specified with their channel id and the distance along that
    channel id.

    Args:
        channel_id (ndarray of int): The channel ids for which to compute weights. Each
            id may occur multiple times, if there are multiple points on that channel
            to compute the weights on. This array must be in ascending order.
        ds (ndarray of float): The location of the points to compute the weights for,
            measured as a distance along the channel. Array must be the same size as
            channel_id and ordered (per channel) in ascending order.
        cs (CrossSectionLocations)
        channels (Channels): Used to lookup the channel geometry

    Returns:
        tuple of cross1, cross2, cross_weight

    See also:
        Channels.get_lines: computes the lines in the correct order
    """
    if not np.in1d(channels.id, cs.channel_id).all():
        missing = set(channels.id) - set(cs.channel_id)
        raise ValueError(f"Channels {missing} have no cross section location set.")

    if len(channel_id) == 0:
        return (
            np.empty((0,), dtype=int),
            np.empty((0,), dtype=int),
            np.empty((0,), dtype=float),
        )

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
    _cumulative = np.cumsum(ds)
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
    line_ch_id = channel_id

    # Fix situations where cs_idx_1 is incorrect
    extrapolate = (cs_ch_id[cs_idx_1] != line_ch_id) | out_of_bounds_1
    cs_idx_1[extrapolate] = cs_idx_2[extrapolate]
    cs_idx_2[extrapolate] = np.clip(cs_idx_2[extrapolate] + 1, None, len(cs) - 1)
    equalize = extrapolate & (cs_ch_id[cs_idx_2] != line_ch_id)
    cs_idx_2[equalize] -= 1

    # Fix situations where cs_idx_2 is incorrect
    extrapolate = (cs_ch_id[cs_idx_2] != line_ch_id) | out_of_bounds_2
    cs_idx_2[extrapolate] = cs_idx_1[extrapolate]
    cs_idx_1[extrapolate] = np.clip(cs_idx_2[extrapolate] - 1, 0, None)
    equalize = extrapolate & (cs_ch_id[cs_idx_1] != line_ch_id)
    cs_idx_1[equalize] += 1

    # Map index to id and create the array that matches the input channel_id and ds
    cross1 = cs.id[cs_sorter][cs_idx_1]
    cross2 = cs.id[cs_sorter][cs_idx_2]

    # Compute the weights. For each line, we have 3 times ds:
    # 1. the ds of the CrossSectionLocation before it (ds_1)
    # 2. the ds of the CrossSectionLocation after it (ds_2)
    # 3. the ds of the midpoint (velocity point) (midpoint_ds)
    #
    # The weight is calculated such that a value at midpoint can be interpolated
    # as:   weight * v_1 + (1 - weight) * v_2
    ds_1 = cs_ds[cs_idx_1]
    ds_2 = cs_ds[cs_idx_2]
    with np.errstate(divide="ignore", invalid="ignore"):
        # this transforms 1 / 0 to inf and 0 / 0 to nan without warning
        cross_weight = (ds_2 - midpoint_ds) / (ds_2 - ds_1)
    cross_weight[~np.isfinite(cross_weight)] = 1.0

    return cross1, cross2, cross_weight


def interpolate(cross1, cross2, cross_weight, cs, cs_attr):
    """Interpolate an attribute of CrossSectionLocation using known weights.

    Args:
        cross1 (ndarray of int): see compute_weights
        cross2 (ndarray of int): see compute_weights
        cross_weight (ndarray of float): see compute_weights
        cs (CrossSectionLocations)
        cs_attr ({"reference_level", "bank_level"}): see compute_weights

    Returns:
        an array of the same shape cross1 containing the interpolated values

    See also:
        compute_weights: computes the interpolation weights
    """
    values = getattr(cs, cs_attr)
    left = values.take(cs.id_to_index(cross1))
    right = values.take(cs.id_to_index(cross2))
    return cross_weight * left + (1 - cross_weight) * right


def compute_bottom_level(channel_id, ds, cs, channels):
    """Compute the bottom level by interpolating/extrapolating between cross sections

    This can be used at nodes (for dmax) or at line centres (for dpumax).

    Args:
        channel_id (ndarray of int): see compute_weights
        ds (ndarray of float): see compute_weights
        cs (CrossSectionLocations): the reference_level is inter/extrapolated
        channels (Channels): see compute_weights

    Returns:
        an array of the same shape as channel_id containing the interpolated values

    See also:
        compute_weights: computes the interpolation weights
    """
    cross1, cross2, cross_weight = compute_weights(channel_id, ds, cs, channels)
    return interpolate(cross1, cross2, cross_weight, cs, "reference_level")


def fix_dpumax(lines, nodes, cs, allow_nan=False):
    """Fix the line bottom levels (dpumax) for channels that have no added nodes.

    The new value is the reference_level of the channel's cross section location
    *only* if that is higher than the already present dpumax. If the channel has
    multiple cs locations, inter/extrapolate on the line midpoint.

    This should be called *after* assigning dpumax to all lines and *after* computing
    the cross section weights.

    Args:
        nodes (Nodes)
        lines (Lines): the dpumax is adjusted where necessary. needs the cross1, cross2
            and cross_weight attributes (see compute_weights)
        cs (CrossSectionLocations): the reference_level is taken
    """
    # find the channel lines that connect 2 connection nodes
    line_idx = np.where(lines.content_type == ContentType.TYPE_V2_CHANNEL)[0]
    node_idx = nodes.id_to_index(lines.line[line_idx])
    is_cn = nodes.content_type[node_idx] == ContentType.TYPE_V2_CONNECTION_NODES
    line_idx = line_idx[is_cn.all(axis=1)]

    # collect the associated crosssection reference_levels and interpolate
    left = cs.reference_level.take(cs.id_to_index(lines.cross1[line_idx]))
    right = cs.reference_level.take(cs.id_to_index(lines.cross2[line_idx]))
    weights = lines.cross_weight[line_idx]
    new_dpumax = weights * left + (1 - weights) * right

    # set the new dpumax, including only the lines that have a lower dpumax
    mask = lines.dpumax[line_idx] < new_dpumax
    lines.dpumax[line_idx[mask]] = new_dpumax[mask]
