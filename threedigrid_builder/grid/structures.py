from threedigrid_builder.base import array_of
from threedigrid_builder.base import Lines
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import FrictionType
from threedigrid_builder.grid import linear

import itertools
import numpy as np
import pygeos


__all__ = ["Culverts", "Orifices", "Weirs"]


class Culvert:  # NL: duiker
    id: int
    code: str
    the_geom: pygeos.Geometry
    dist_calc_points: float
    connection_node_start_id: int
    connection_node_end_id: int
    calculation_type: CalculationType
    cross_section_definition_id: int
    invert_level_start_point: float
    invert_level_end_point: float
    discharge_coefficient_negative: float
    discharge_coefficient_positive: float
    friction_type: FrictionType
    friction_value: float
    # zoom_category
    # display_name


@array_of(Culvert)
class Culverts(linear.BaseLinear):
    def set_geometries(self, connection_nodes):
        """Compute culvert geometries from the connection nodes where necessary"""
        has_no_geom = pygeos.is_missing(self.the_geom)
        if not has_no_geom.any():
            return

        # construct the culvert geometries
        points_1 = connection_nodes.the_geom[
            connection_nodes.id_to_index(self.connection_node_start_id[has_no_geom])
        ]
        points_2 = connection_nodes.the_geom[
            connection_nodes.id_to_index(self.connection_node_end_id[has_no_geom])
        ]
        coordinates = np.empty((np.count_nonzero(has_no_geom), 2, 2))
        coordinates[:, 0, 0] = pygeos.get_x(points_1)
        coordinates[:, 0, 1] = pygeos.get_y(points_1)
        coordinates[:, 1, 0] = pygeos.get_x(points_2)
        coordinates[:, 1, 1] = pygeos.get_y(points_2)
        self.the_geom[has_no_geom] = pygeos.linestrings(coordinates)

    def interpolate_nodes(self, *args, **kwargs):
        """Compute interpolated nodes for culverts.

        See also:
            BaseLinear.interpolate_nodes
        """
        if pygeos.is_missing(self.the_geom).any():
            raise ValueError(
                "Culverts found without a geometry. Call set_geometries first."
            )
        nodes = super().interpolate_nodes(*args, **kwargs)
        nodes.content_type[:] = ContentType.TYPE_V2_CULVERT
        return nodes

    def get_lines(self, *args, **kwargs):
        """Compute the grid lines for the culverts.

        See also:
            BaseLinear.get_lines
        """
        lines = super().get_lines(*args, **kwargs)
        lines.content_type[:] = ContentType.TYPE_V2_CULVERT
        return lines


class _WeirOrifice:  # NL: stuw / doorlaat
    id: int
    code: str
    connection_node_start_id: int
    connection_node_end_id: int
    crest_level: float
    crest_type: CalculationType
    cross_section_definition_id: int
    discharge_coefficient_negative: float
    discharge_coefficient_positive: float
    friction_type: FrictionType
    friction_value: float
    # zoom_category
    # display_name
    # sewerage


@array_of(_WeirOrifice)
class _WeirOrifices:
    content_type = None  # subclasses need to define this

    def get_lines(self, connection_nodes, line_id_counter, connection_node_offset=0):
        """Convert weirs (or orifices) into lines

        Args:
            connection_nodes (ConnectionNodes): used to map ids to indices
            line_id_counter (iterable): an iterable yielding integers
            connection_node_offset (int): offset to give connection node
              indices in the returned lines.line. Default 0.

        Returns:
            Lines with data in the following columns:
            - id: counter generated from line_id_counter
            - line: 2 node ids per line
            - content_type: the content type of the structure (weir or orifice)
            - content_pk: the id of the structure from which this line originates
            - kcu: the crest_type of the structure
            - dpumax: the crest_level of the structure
        """
        # map connection node IDs to node indices
        cn_start_idx = connection_nodes.id_to_index(self.connection_node_start_id)
        cn_end_idx = connection_nodes.id_to_index(self.connection_node_end_id)
        line = np.array([cn_start_idx, cn_end_idx]) + connection_node_offset
        return Lines(
            id=itertools.islice(line_id_counter, len(self)),
            line=line,
            content_type=self.content_type,
            content_pk=self.id,
            kcu=self.crest_type,  # implicitly converts CalculationType -> LineType
            dpumax=self.crest_level,
        )


class Orifices(_WeirOrifices):
    content_type = ContentType.TYPE_V2_ORIFICE


class Weirs(_WeirOrifices):
    content_type = ContentType.TYPE_V2_WEIR
