import itertools

import numpy as np
import shapely

from threedigrid_builder.base import Array, Lines
from threedigrid_builder.constants import CalculationType, ContentType, FrictionType
from threedigrid_builder.grid import linear

__all__ = ["Culverts", "Orifices", "Weirs"]


class Culvert:  # NL: duiker
    id: int
    code: str
    the_geom: shapely.Geometry
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
    zoom_category: int
    display_name: str


class Culverts(Array[Culvert], linear.BaseLinear):
    content_type = ContentType.TYPE_V2_CULVERT

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
    zoom_category: int
    display_name: str
    sewerage: int


class WeirOrifices(Array[WeirOrifice]):
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
            - cross_id1 & cross_id2: the id of the cross section definition
            - cross_weight: 1.0 (which means that cross_id2 should be ignored)
            - frict_type1 & frict_type2: the friction type (both are equal)
            - frict_value1 & frict_value2: the friction value (both are equal)
            - invert_level_start_point: the crest_level of the structure
            - invert_level_end_point: the crest_level of the structure
            - discharge_coefficient_positive and _negative: taken from the structure
            - display_name: user defined label for display purposes
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
            cross_id1=self.cross_section_definition_id,
            cross_id2=self.cross_section_definition_id,
            cross_weight=1.0,
            frict_type1=self.friction_type,
            frict_value1=self.friction_value,
            frict_type2=self.friction_type,
            frict_value2=self.friction_value,
            invert_level_start_point=self.crest_level,
            invert_level_end_point=self.crest_level,
            discharge_coefficient_positive=self.discharge_coefficient_positive,
            discharge_coefficient_negative=self.discharge_coefficient_negative,
            display_name=self.display_name,
            connection_node_start_id=self.connection_node_start_id,
            connection_node_end_id=self.connection_node_end_id,
            crest_level=self.crest_level,
            crest_type=self.crest_type,
            sewerage=self.sewerage,
        )


class Orifices(WeirOrifices):
    content_type = ContentType.TYPE_V2_ORIFICE


class Weirs(WeirOrifices):
    content_type = ContentType.TYPE_V2_WEIR
