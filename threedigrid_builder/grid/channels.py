import numpy as np
import shapely

from threedigrid_builder.base import Array
from threedigrid_builder.constants import CalculationType, ContentType
from threedigrid_builder.grid import linear
from threedigrid_builder.grid.cross_section_locations import compute_bottom_level

__all__ = ["Channels"]


class Channel:
    id: int
    code: str
    the_geom: shapely.Geometry
    dist_calc_points: float
    connection_node_start_id: int
    connection_node_end_id: int
    calculation_type: CalculationType
    display_name: str
    zoom_category: int
    exchange_thickness: float
    hydraulic_conductivity_out: float
    hydraulic_conductivity_in: float


class Channels(Array[Channel], linear.BaseLinear):
    content_type = ContentType.TYPE_V2_CHANNEL

    def is_closed(self, content_pk):
        """Whether objects are 'closed' or 'open water'.

        This is relevant for 1D-2D connections.
        """
        return np.full(len(content_pk), False, dtype=bool)

    def get_1d2d_exchange_levels(self, content_pk, s1d, locations, **kwargs):
        """Compute the exchange level (dpumax) for 1D-2D flowlines.

        For channels 1D-2D exchange levels are interpolated between crosssection
        location bank levels.

        Args:
            content_pk (array of int): object ids for which to compute levels
            s1d (array of int): distance on the objects where to compute levels
            locations (CrossSectionLocations): for the bank_levels

        Returns:
            exchange levels a.k.a. dpumax (array of float)
        """
        return compute_bottom_level(content_pk, s1d, locations, self, "bank_level")

    @property
    def has_groundwater_exchange(self):
        with np.errstate(invalid="ignore"):
            return (
                (self.exchange_thickness > 0)
                & np.isfinite(self.hydraulic_conductivity_out)
                & np.isfinite(self.hydraulic_conductivity_in)
            )
