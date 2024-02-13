import numpy as np
import shapely

from threedigrid_builder.base import Array, PointsOnLine
from threedigrid_builder.constants import FrictionType

__all__ = ["CrossSectionLocations"]


class CrossSectionLocation:
    id: int
    code: str
    the_geom: shapely.Geometry
    definition_id: id  # refers to CrossSectionDefinition
    channel_id: id  # refers to Channel
    reference_level: float
    bank_level: float
    friction_type: FrictionType
    friction_value: float
    vegetation_stem_density: float
    vegetation_stem_diameter: float
    vegetation_height: float
    vegetation_drag_coefficient: float


class CrossSectionLocations(Array[CrossSectionLocation]):
    def apply_to_lines(self, lines, channels, extrapolate=False):
        """Apply cross section locations to lines

        Args:
            lines (Lines): Changed inplace
            channels (Channels)
            extrapolate (bool): Whether to allow extrapolation. Extrapolation may occur
                when there are multiple cross section locations and some line centers
                are not in between them. When turned off, the values of those lines
                become equal to the values of the closest cross section location.
                Default False.

        The following lines attributes are changed inplace:

        - cross_id1: the id of the first cross section definition
        - cross_id2: the id of the second cross section definition
        - cross_weight: the weight of the first cross section definition
        - frict_type1: the friction type of the first cross section location
        - frict_type2: the friction type of the second cross section location
        - frict_value1: the friction value of the first cross section location
        - frict_value2: the friction value of the second cross section location
        - veg_coef1: the product of vegetation properties of the first cross section location (except height)
        - veg_coef2: the product of vegetation properties of the second cross section location (except height)
        - veg_height1: the vegetation height of the first cross section
        - veg_height2: the vegetation height of the second cross section
        - invert_level_start_point: 'reference_level' interpolated at the line end
        - invert_level_end_point: 'reference_level' interpolated at the line start
        - dpumax: the largest of the two invert levels
        """
        # Mask the lines to only the Channel lines
        cross_loc1, cross_loc2, cross_weight = compute_weights(
            lines.content_pk,
            lines.s1d,
            self,
            channels,
            extrapolate=extrapolate,
        )

        idx1 = self.id_to_index(cross_loc1)
        idx2 = self.id_to_index(cross_loc2)
        lines.cross_id1 = self.definition_id[idx1]
        lines.cross_id2 = self.definition_id[idx2]
        lines.cross_weight = cross_weight

        lines.frict_type1 = self.friction_type[idx1]
        lines.frict_type2 = self.friction_type[idx2]
        lines.frict_value1 = self.friction_value[idx1]
        lines.frict_value2 = self.friction_value[idx2]

        lines.veg_coef1 = (
            self.vegetation_stem_density[idx1]
            * self.vegetation_stem_diameter[idx1]
            * self.vegetation_drag_coefficient[idx1]
        )
        lines.veg_coef1[np.isnan(lines.veg_coef1)] = 0.0

        lines.veg_coef2 = (
            self.vegetation_stem_density[idx2]
            * self.vegetation_stem_diameter[idx2]
            * self.vegetation_drag_coefficient[idx2]
        )
        lines.veg_coef2[np.isnan(lines.veg_coef2)] = 0.0
        lines.veg_height1 = self.vegetation_height[idx1]
        lines.veg_height1[np.isnan(lines.veg_height1)] = 0.0
        lines.veg_height2 = self.vegetation_height[idx2]
        lines.veg_height2[np.isnan(lines.veg_height2)] = 0.0

        # Compute invert levels and start and end
        lines.invert_level_start_point = compute_bottom_level(
            lines.content_pk,
            lines.s1d - (lines.ds1d / 2),
            self,
            channels,
            extrapolate=extrapolate,
        )
        lines.invert_level_end_point = compute_bottom_level(
            lines.content_pk,
            lines.s1d + (lines.ds1d / 2),
            self,
            channels,
            extrapolate=extrapolate,
        )
        lines.dpumax = np.maximum(
            lines.invert_level_start_point, lines.invert_level_end_point
        )


def compute_weights(channel_id, points, cs, channels, extrapolate):
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
        extrapolate (bool): Whether to allow the weight to go outside of [0, 1]

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

    cs_points = PointsOnLine.from_geometries(
        channels.linestrings,
        cs.the_geom,
        channels.id_to_index(cs.channel_id),
        content_pk=cs.id,
    )
    queried_points = PointsOnLine.from_s1d(
        channels.linestrings, points, channels.id_to_index(channel_id)
    )

    # Find what CS location comes after and before each midpoint
    cs_idx_1, cs_idx_2 = cs_points.neighbours(queried_points)

    # Map index to id and create the array that matches the input channel_id and ds
    cross_loc1 = cs_points.content_pk[cs_idx_1]
    cross_loc2 = cs_points.content_pk[cs_idx_2]

    # Compute the weights. For each line, we have 3 times s:
    # 1. the s of the CrossSectionLocation before it (cs_idx_1)
    # 2. the s of the CrossSectionLocation after it (cs_idx_2)
    # 3. the s of the midpoint (velocity point) (queried_points)
    #
    # The weight is calculated such that a value at midpoint can be interpolated
    # as:   weight * v_1 + (1 - weight) * v_2
    s_1 = cs_points.s1d[cs_idx_1]
    s_2 = cs_points.s1d[cs_idx_2]
    with np.errstate(divide="ignore", invalid="ignore"):
        # this transforms 1 / 0 to inf and 0 / 0 to nan without warning
        cross_weight = (s_2 - queried_points.s1d) / (s_2 - s_1)
    cross_weight[~np.isfinite(cross_weight)] = 1.0

    if not extrapolate:
        cross_weight = np.clip(cross_weight, 0.0, 1.0)

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


def compute_bottom_level(
    channel_id, ds, cs, channels, cs_attr="reference_level", extrapolate=False
):
    """Compute levels by interpolating/extrapolating between cross sections

    This can be used at nodes (for dmax) or at line centres (for dpumax).

    Args:
        channel_id (ndarray of int): see compute_weights
        ds (ndarray of float): see compute_weights
        cs (CrossSectionLocations): the reference_level is inter/extrapolated
        channels (Channels): see compute_weights
        cs_attr ({"reference_level", "bank_level"})
        extrapolate (bool): Whether to allow extrapolation. Extrapolation may occur
            when there are multiple cross section locations and some points (ds)
            are not in between them. When turned off, the values at those points
            become equal to the values of the closest cross section locations.
            Default False.

    Returns:
        an array of the same shape as channel_id containing the interpolated values

    See also:
        compute_weights: computes the interpolation weights
    """
    cross_loc1, cross_loc2, cross_weight = compute_weights(
        channel_id, ds, cs, channels, extrapolate=extrapolate
    )
    return interpolate(cross_loc1, cross_loc2, cross_weight, cs, cs_attr)
