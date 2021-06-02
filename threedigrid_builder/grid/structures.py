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
    content_type = ContentType.TYPE_V2_CULVERT


class WeirOrifice:  # NL: stuw / doorlaat
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


@array_of(WeirOrifice)
class WeirOrifices:
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
        line = np.empty((len(self), 2), dtype=np.int32, order="F")
        line[:, 0] = connection_nodes.id_to_index(self.connection_node_start_id)
        line[:, 1] = connection_nodes.id_to_index(self.connection_node_end_id)
        line += connection_node_offset
        return Lines(
            id=itertools.islice(line_id_counter, len(self)),
            line=line,
            content_type=self.content_type,
            content_pk=self.id,
            kcu=self.crest_type,  # implicitly converts CalculationType -> LineType
            dpumax=self.crest_level,
        )


class Orifices(WeirOrifices):
    content_type = ContentType.TYPE_V2_ORIFICE


class Weirs(WeirOrifices):
    content_type = ContentType.TYPE_V2_WEIR
