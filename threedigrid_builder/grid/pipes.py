import numpy as np
import shapely

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
    the_geom: shapely.Geometry
    zoom_category: int
    display_name: str
    # profile_num
    # original_length
    material: int
    exchange_thickness: float
    hydraulic_conductivity_out: float
    hydraulic_conductivity_in: float


class Pipes(Array[Pipe], linear.BaseLinear):
    content_type = ContentType.TYPE_V2_PIPE

    def is_closed(self, content_pk):
        """Whether objects are 'closed' or 'open water'.

        This is relevant for 1D-2D connections.
        """
        return np.full(len(content_pk), True, dtype=bool)

    def get_1d2d_exchange_levels(self, content_pk, s1d, connection_nodes, **kwargs):
        """Compute the exchange level (dpumax) for 1D-2D flowlines.

        For pipes and culverts 1D2D exchange levels are interpololated between
        connection node drain levels.

        Args:
            content_pk (array of int): object ids for which to compute levels
            s1d (array of int): distance on the objects where to compute levels
            connection_nodes (ConnectionNodes): for the drain_levels

        Returns:
            exchange levels a.k.a. dpumax (array of float)
        """
        return self.compute_drain_level(
            ids=content_pk,
            s=s1d,
            connection_nodes=connection_nodes,
        )

    @property
    def has_groundwater_exchange(self):
        with np.errstate(invalid="ignore"):
            return (
                (self.exchange_thickness > 0)
                & np.isfinite(self.hydraulic_conductivity_out)
                & np.isfinite(self.hydraulic_conductivity_in)
            )
