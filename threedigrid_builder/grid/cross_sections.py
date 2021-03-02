from threedigrid_builder.base import array_of
from threedigrid_builder.constants import CrossSectionShape
from threedigrid_builder.constants import FrictionType

import pygeos


__all__ = ["CrossSectionLocations", "CrossSectionDefinitions"]


class CrossSectionLocation:
    id: int
    # code: str  unused?
    the_geom: pygeos.Geometry
    definition_id: id  # refers to CrossSectionDefinition
    channel_id: id  # refers to Channel
    reference_level: float
    bank_level: float
    friction_type: FrictionType
    friction_value: float


@array_of(CrossSectionLocation)
class CrossSectionLocations:
    pass


class CrossSectionDefinition:
    id: int
    # code: str  unused?
    shape: CrossSectionShape
    height: float
    width: float


@array_of(CrossSectionDefinition)
class CrossSectionDefinitions:
    pass
