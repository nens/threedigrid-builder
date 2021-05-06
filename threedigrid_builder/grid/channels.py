from threedigrid_builder.base import array_of
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.grid import linear

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
    def interpolate_nodes(self, *args, **kwargs):
        """Compute interpolated nodes for channels.

        See also:
            BaseLinear.interpolate_nodes
        """
        nodes = super().interpolate_nodes(*args, **kwargs)
        nodes.content_type[:] = ContentType.TYPE_V2_CHANNEL
        return nodes

    def get_lines(self, *args, **kwargs):
        """Compute the grid lines for the channels.

        See also:
            BaseLinear.get_lines
        """
        lines = super().get_lines(*args, **kwargs)
        lines.content_type[:] = ContentType.TYPE_V2_CHANNEL
        return lines
