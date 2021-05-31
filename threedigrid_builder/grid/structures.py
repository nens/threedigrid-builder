from threedigrid_builder.base import array_of
from threedigrid_builder.constants import CalculationType
from threedigrid_builder.constants import FrictionType

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
class Culverts:
    pass


class Orifice:
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
    the_geom: pygeos.Geometry
    # zoom_category
    # display_name
    # sewerage


@array_of(Orifice)
class Orifices:
    pass


class Weir:  # NL: stuw
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
    the_geom: pygeos.Geometry
    # zoom_category
    # display_name
    # sewerage


@array_of(Weir)
class Weirs:
    pass
