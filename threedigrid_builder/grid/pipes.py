from threedigrid_builder.base import array_of
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import FrictionType
from threedigrid_builder.constants import NodeType
from threedigrid_builder.constants import SewerageType
from .geo_utils import segmentize

import itertools
import numpy as np
import pygeos


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
    # zoom_category
    # display_name
    # profile_num
    # original_length
    # material


@array_of(Pipe)
class Pipes:
    def interpolate_nodes(
        self, node_id_counter, global_dist_calc_points, connection_nodes
    ):
        """Compute interpolated nodes for pipes.

        The distance between the nodes is set by pipes.dist_calc_points. If this is not
        available, the global setting is used. The pipe geometry is retrieved from its
        two connection nodes.

        Args:
            node_id_counter (iterable): an iterable yielding integers
            global_dist_calc_points (float): Default node interdistance.
            connection_nodes (ConnectionNodes)

        Returns:
            tuple of nodes (Nodes), segment_size (ndarray)
            nodes has data in the following columns:
            - id: 0-based counter generated here
            - coordinates
            - content_type: ContentType.TYPE_V2_PIPE
            - content_pk: the id of the Pipe from which this node originates
            - node_type: NodeType.NODE_1D_NO_STORAGE
        """
        # load data
        dists = self.dist_calc_points.copy()  # copy because of inplace edits

        # insert default dist_calc_points where necessary
        dists[~np.isfinite(dists)] = global_dist_calc_points
        dists[dists <= 0] = global_dist_calc_points

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
        geometries = pygeos.linestrings(coordinates)

        # compute number of nodes to add per pipe
        _, points, _, index, segment_size = segmentize(geometries, dists)

        nodes = Nodes(
            id=itertools.islice(node_id_counter, len(points)),
            coordinates=pygeos.get_coordinates(points),
            content_type=ContentType.TYPE_V2_PIPE,
            content_pk=self.index_to_id(index),
            node_type=NodeType.NODE_1D_NO_STORAGE,
            calculation_type=self.calculation_type[index],
        )
        return nodes, segment_size
