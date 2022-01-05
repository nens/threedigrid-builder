from threedigrid_builder.base import array_of
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import linear
from threedigrid_builder.grid.cross_section_locations import (
    compute_bottom_level,
)

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
    display_name: str
    zoom_category: int


@array_of(Channel)
class Channels(linear.BaseLinear):
    content_type = ContentType.TYPE_V2_CHANNEL

    def get_1d2d_properties(self, nodes, node_idx, locations):
        """Compute properties (is_closed, dpumax) of 1D-2D channel flowlines.

        Args:
            nodes (Nodes): All nodes
            node_idx (array of int): indices into nodes for which to compute properties
            locations (CrossSectionLocations): for the bank_levels

        Returns:
            tuple of:
            - is_closed (bool): always False
            - dpumax (array of float): interpolated between CS location bank_levels
        """
        # dpumax is interpolated between cross section location bank levels
        dpumax = compute_bottom_level(
            nodes.content_pk[node_idx],
            nodes.s1d[node_idx],
            locations,
            self,
            "bank_level",
        )

        return False, dpumax
