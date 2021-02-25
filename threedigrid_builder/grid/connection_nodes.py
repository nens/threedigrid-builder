from threedigrid_builder.base import array_of
from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from threedigrid_builder.constants import ContentType

import numpy as np
import pygeos


__all__ = ["ConnectionNodes"]


class ConnectionNode:
    id: int
    the_geom: pygeos.Geometry
    code: str
    storage_area: float


@array_of(ConnectionNode)
class ConnectionNodes:
    pass
