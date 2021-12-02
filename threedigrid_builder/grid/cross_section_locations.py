from threedigrid_builder.base import array_of
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import FrictionType
from threedigrid_builder.exceptions import SchematisationError

import numpy as np
import pygeos


__all__ = ["CrossSectionLocations"]


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
    def apply_to_lines(self, lines, channels, definitions):
        """Apply cross section locations to lines

        The following lines attributes are changed inplace:

        - cross1: the id of the first cross section definition
        - cross2: the id of the second cross section definition
        - cross_weight: the weight of the first cross section definition
        - invert_level_start_point: 'reference_level' interpolated at the line end
        - invert_level_end_point: 'reference_level' interpolated at the line start
        """
        # Mask the lines to only the Channel lines
        cross_loc1, cross_loc2, cross_weight = compute_weights(
            lines.content_pk,
            lines.s1d,
            self,
            channels,
        )

        idx1 = self.id_to_index(cross_loc1)
        idx2 = self.id_to_index(cross_loc2)
        # Fill cross1 and cross2 by mapping to CrossSectionDefinitions
        try:
            lines.cross1 = definitions.id_to_index(
                self.definition_id[idx1],
                check_exists=True,
            )
            lines.cross2 = definitions.id_to_index(
                self.definition_id[idx2],
                check_exists=True,
            )
        except KeyError as e:
            raise SchematisationError(
                f"Cross section locations {self.id[e.indices].tolist()} refer to "
                f"non-existing cross section definitions."
            )
        lines.cross_weight = cross_weight

        lines.frict_type1 = self.friction_type[idx1]
        lines.frict_type2 = self.friction_type[idx2]
        lines.frict_value1 = self.friction_value[idx1]
        lines.frict_value2 = self.friction_value[idx2]

        # Compute invert levels and start and end
        lines.invert_level_start_point = compute_bottom_level(
            lines.content_pk,
            lines.s1d - (lines.ds1d / 2),
            self,
            channels,
        )
        lines.invert_level_end_point = compute_bottom_level(
            lines.content_pk,
            lines.s1d + (lines.ds1d / 2),
            self,
            channels,
        )


def compute_weights(channel_id, points, cs, channels):
    """Compute cross section weights for points on channels.

    Points on channels are specified with their channel id and the distance along that
    channel id.

    Args:
        channel_id (ndarray of int): The channel ids for which to compute weights. Each
            id may occur multiple times, if there are multiple points on that channel
            to compute the weights on. This array must be in ascending order.
        points (ndarray of float): The location of the points to compute the weights for,
            measured as a distance along the channel. Array must be the same size as
            channel_id and ordered (per channel) in ascending order.
        cs (CrossSectionLocations)
        channels (Channels): Used to lookup the channel geometry

    Returns:
        tuple of cross_loc1, cross_loc2, cross_weight

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

    if np.any(~np.isfinite(points)):
        raise ValueError("NaN values encountered in points")

    cs_channel_idx = channels.id_to_index(cs.channel_id)
    points_channel_idx = channels.id_to_index(channel_id)

    # Locate the position measured along the channel (ds) of the CS locations
    cs_s = pygeos.line_locate_point(channels.the_geom[cs_channel_idx], cs.the_geom)

    # To each cs_s, add the length of all channels before it
    # the lengths are increased by a small number to mitigate overlapping start/ends
    ch_cum_length = np.cumsum(pygeos.length(channels.the_geom) + 1e-6)
    ch_cum_length = np.roll(ch_cum_length, 1)
    ch_cum_length[0] = 0.0
    cs_s += ch_cum_length[cs_channel_idx]

    # Make cs_s monotonically increasing and update cs_sorter to match it
    cs_sorter = np.argsort(cs_s)
    cs_s = cs_s[cs_sorter]
    cs_channel_idx = cs_channel_idx[cs_sorter]

    # Compute the ds of the line points, cumulative over all channels
    points_cum = points + ch_cum_length[points_channel_idx]

    # Find what CS location comes after and before each midpoint
    cs_idx_2 = np.searchsorted(cs_s, points_cum)
    cs_idx_1 = cs_idx_2 - 1
    out_of_bounds_1 = cs_idx_1 < 0
    out_of_bounds_2 = cs_idx_2 >= len(cs)
    cs_idx_1[out_of_bounds_1] = 0
    cs_idx_2[out_of_bounds_2] = len(cs) - 1

    # Extrapolation: if a point does not have a cross section location
    # to its left (cs_idx_1 is invalid), assign to it the 2 cross_idx
    # to its right. And vice versa for cs_idx_2 that is invalid. Note
    # that because every channel has at least 1 crosssection, one the two
    # must always be valid. This is not checked here.
    # Equalization: if a newly assigned crosssection location would not
    # belong to the correct channel, we have only 1 crossection location
    # in the channel so we assign the same cross_idx to both 1 and 2.

    # Fix situations where cs_idx_1 is incorrect
    extrapolate = (cs_channel_idx[cs_idx_1] != points_channel_idx) | out_of_bounds_1
    cs_idx_1[extrapolate] = cs_idx_2[extrapolate]
    cs_idx_2[extrapolate] = np.clip(cs_idx_2[extrapolate] + 1, None, len(cs) - 1)
    equalize = extrapolate & (cs_channel_idx[cs_idx_2] != points_channel_idx)
    cs_idx_2[equalize] -= 1

    # Fix situations where cs_idx_2 is incorrect
    extrapolate = (cs_channel_idx[cs_idx_2] != points_channel_idx) | out_of_bounds_2
    cs_idx_2[extrapolate] = cs_idx_1[extrapolate]
    cs_idx_1[extrapolate] = np.clip(cs_idx_2[extrapolate] - 1, 0, None)
    equalize = extrapolate & (cs_channel_idx[cs_idx_1] != points_channel_idx)
    cs_idx_1[equalize] += 1

    # Map index to id and create the array that matches the input channel_id and ds
    cross_loc1 = cs.id[cs_sorter][cs_idx_1]
    cross_loc2 = cs.id[cs_sorter][cs_idx_2]

    # Compute the weights. For each line, we have 3 times s:
    # 1. the s of the CrossSectionLocation before it (s_1)
    # 2. the s of the CrossSectionLocation after it (s_2)
    # 3. the s of the midpoint (velocity point) (points)
    #
    # The weight is calculated such that a value at midpoint can be interpolated
    # as:   weight * v_1 + (1 - weight) * v_2
    s_1 = cs_s[cs_idx_1]
    s_2 = cs_s[cs_idx_2]
    with np.errstate(divide="ignore", invalid="ignore"):
        # this transforms 1 / 0 to inf and 0 / 0 to nan without warning
        cross_weight = (s_2 - points_cum) / (s_2 - s_1)
    cross_weight[~np.isfinite(cross_weight)] = 1.0

    return cross_loc1, cross_loc2, cross_weight


def interpolate(cross_loc1, cross_loc2, cross_weight, cs, cs_attr="reference_level"):
    """Interpolate an attribute of CrossSectionLocation using known weights.

    Args:
        cross_loc1 (ndarray of int): see compute_weights
        cross_loc2 (ndarray of int): see compute_weights
        cross_weight (ndarray of float): see compute_weights
        cs (CrossSectionLocations)
        cs_attr ({"reference_level", "bank_level"})

    Returns:
        an array of the same shape cross_loc1 containing the interpolated values

    See also:
        compute_weights: computes the interpolation weights
    """
    values = getattr(cs, cs_attr)
    left = values.take(cs.id_to_index(cross_loc1))
    right = values.take(cs.id_to_index(cross_loc2))
    return cross_weight * left + (1 - cross_weight) * right


def compute_bottom_level(channel_id, ds, cs, channels, cs_attr="reference_level"):
    """Compute levels by interpolating/extrapolating between cross sections

    This can be used at nodes (for dmax) or at line centres (for dpumax).

    Args:
        channel_id (ndarray of int): see compute_weights
        ds (ndarray of float): see compute_weights
        cs (CrossSectionLocations): the reference_level is inter/extrapolated
        channels (Channels): see compute_weights
        cs_attr ({"reference_level", "bank_level"})

    Returns:
        an array of the same shape as channel_id containing the interpolated values

    See also:
        compute_weights: computes the interpolation weights
    """
    cross_loc1, cross_loc2, cross_weight = compute_weights(channel_id, ds, cs, channels)
    return interpolate(cross_loc1, cross_loc2, cross_weight, cs, cs_attr)


def fix_dpumax(lines, nodes):
    """Fix the line bottom levels (dpumax) for channels that have no added nodes.

    The new value is the reference_level of the channel's cross section location
    *only* if that is higher than the already present dpumax. If the channel has
    multiple cs locations, inter/extrapolate on the line midpoint.

    This should be called *after* assigning dpumax to all lines and *after* computing
    the invert levels of the lines.

    Args:
        nodes (Nodes)
        lines (Lines): the dpumax is adjusted where necessary. Needs the
            invert_level_start_point and invert_level_end_point attributes.
    """
    # find the channel lines that connect 2 connection nodes
    line_idx = np.where(lines.content_type == ContentType.TYPE_V2_CHANNEL)[0]
    node_idx = nodes.id_to_index(lines.line[line_idx])
    is_cn = nodes.content_type[node_idx] == ContentType.TYPE_V2_CONNECTION_NODES
    line_idx = line_idx[is_cn.all(axis=1)]

    # average between the invert levels
    new_dpumax = (
        lines.invert_level_start_point[line_idx]
        + lines.invert_level_end_point[line_idx]
    ) / 2

    # set the new dpumax, including only the lines that have a lower dpumax
    mask = lines.dpumax[line_idx] < new_dpumax
    lines.dpumax[line_idx[mask]] = new_dpumax[mask]
