from threedigrid_builder.base import array_of
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import FrictionType
from threedigrid_builder.constants import SewerageType
from threedigrid_builder.grid import linear

import numpy as np
import pygeos


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
    sewerage_type: SewerageType
    friction_type: FrictionType
    friction_value: float
    the_geom: pygeos.Geometry
    # zoom_category
    # display_name
    # profile_num
    # original_length
    # material


@array_of(Pipe)
class Pipes(linear.BaseLinear):
    def set_geometries(self, connection_nodes):
        """Compute pipe geometries from the connection node locations"""
        # construct the pipe geometries
        points_1 = connection_nodes.the_geom[
            connection_nodes.id_to_index(self.connection_node_start_id)
        ]
        points_2 = connection_nodes.the_geom[
            connection_nodes.id_to_index(self.connection_node_end_id)
        ]
        coordinates = np.empty((len(self), 2, 2))
        coordinates[:, 0, 0] = pygeos.get_x(points_1)
        coordinates[:, 0, 1] = pygeos.get_y(points_1)
        coordinates[:, 1, 0] = pygeos.get_x(points_2)
        coordinates[:, 1, 1] = pygeos.get_y(points_2)
        self.the_geom = pygeos.linestrings(coordinates)

    def interpolate_nodes(self, *args, **kwargs):
        """Compute interpolated nodes for pipes.

        See also:
            BaseLinear.interpolate_nodes
        """
        if pygeos.is_missing(self.the_geom).any():
            raise ValueError(
                "Pipes found without a geometry. Call set_geometries first."
            )
        nodes = super().interpolate_nodes(*args, **kwargs)
        nodes.content_type[:] = ContentType.TYPE_V2_PIPE
        return nodes

    def get_lines(self, *args, **kwargs):
        """Compute the grid lines for the pipes.

        See also:
            BaseLinear.get_lines
        """
        lines = super().get_lines(*args, **kwargs)
        lines.content_type[:] = ContentType.TYPE_V2_PIPE
        return lines


def compute_bottom_level(pipe_id, ds, pipes):
    """Compute the bottom level by interpolating between invert levels

    This function is to be used for interpolated nodes on pipes.

    Args:
        pipe_id (ndarray of int): the id (content_pk) of the pipe
        ds (ndarray of float): the distance from the pipe start to the point
        pipes (Pipes): used to lookup the invert levels

    Returns:
        an array of the same shape as pipe_id and ds containing the interpolated values
    """
    if pygeos.is_missing(pipes.the_geom).any():
        raise ValueError("Pipes found without a geometry. Call set_geometries first.")
    pipe_lengths = pygeos.length(pipes.the_geom)
    pipe_idx = pipes.id_to_index(pipe_id)
    weights = ds / pipe_lengths[pipe_idx]
    if np.any(weights < 0.0) or np.any(weights > 1.0):
        raise ValueError("Encountered pipe nodes outside of the pipe bounds")

    left = pipes.invert_level_start_point[pipe_idx]
    right = pipes.invert_level_end_point[pipe_idx]

    return weights * right + (1 - weights) * left
