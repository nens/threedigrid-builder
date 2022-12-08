import pygeos

from threedigrid_builder.base import Array
from threedigrid_builder.constants import CalculationType, ContentType, FrictionType
from threedigrid_builder.grid import linear

__all__ = ["Pipes"]


class Pipe:
    id: int
    code: str
    dist_calc_points: float
    calculation_type: CalculationType
    connection_node_start_id: int
    connection_node_end_id: int
    cross_section_definition_id: int
    invert_level_start_point: float
    invert_level_end_point: float
    sewerage_type: int
    friction_type: FrictionType
    friction_value: float
    the_geom: pygeos.Geometry
    zoom_category: int
    display_name: str
    # profile_num
    # original_length
    material: int


class Pipes(Array[Pipe], linear.BaseLinear):
    content_type = ContentType.TYPE_V2_PIPE

    def get_1d2d_properties(self, nodes, node_idx, connection_nodes, **kwargs):
        """Compute properties (is_closed, dpumax) of 1D-2D pipe flowlines.

        Args:
            nodes (Nodes): All nodes
            node_idx (array of int): indices into nodes for which to compute properties
            connection_nodes (ConnectionNodes): for the drain_level

        Returns:
            tuple of:
            - is_closed (bool): always True
            - dpumax (array of float): interpolated between CN drain_level
        """
        # dpumax is interpolated between drain levels of adjacent manholes (conn nodes)
        dpumax = self.compute_drain_level(
            ids=nodes.content_pk[node_idx],
            s=nodes.s1d[node_idx],
            connection_nodes=connection_nodes,
        )

        return True, dpumax
