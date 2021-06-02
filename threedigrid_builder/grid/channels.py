from threedigrid_builder.base import array_of
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import linear
from threedigrid_builder.grid.cross_sections import interpolate

import pygeos


__all__ = ["Channels"]


class Channel:
    id: int
    code: str
    the_geom: pygeos.Geometry
    dist_calc_points: float
    connection_node_start_id: int
    connection_node_end_id: int
    calculation_type: CalculationType


@array_of(Channel)
class Channels(linear.BaseLinear):
    content_type = ContentType.TYPE_V2_CHANNEL

    def get_1d2d_properties(self, nodes, node_idx, locations):
        """Compute properties (is_sewerage, dpumax) of 1D-2D flowlines.

        Args:
            nodes (Nodes): All nodes
            node_idx (array of int): indices into nodes for which to compute properties
            locations (CrossSectionLocations): for the bank_levels

        Returns:
            tuple of:
            - is_sewerage (bool): always False
            - dpumax (array of float): interpolated between CS location bank_levels
        """
        # dpumax is interpolated between cross section location bank levels
        dpumax = interpolate(
            nodes.cross1[node_idx],
            nodes.cross2[node_idx],
            nodes.cross_weight[node_idx],
            locations,
            "bank_level",
        )

        return False, dpumax
