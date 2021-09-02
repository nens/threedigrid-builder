from threedigrid_builder.base import array_of
from threedigrid_builder.constants import BoundaryType

import pygeos


__all__ = ["BoundaryConditions1D", "BoundaryConditions2D"]


class BoundaryCondition1D:
    id: int
    boundary_type: BoundaryType
    connection_node_id: int


@array_of(BoundaryCondition1D)
class BoundaryConditions1D:
    pass


class BoundaryCondition2D:
    id: int
    boundary_type: BoundaryType
    the_geom: pygeos.Geometry


@array_of(BoundaryCondition2D)
class BoundaryConditions2D:
    pass
