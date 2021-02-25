from threedigrid_builder.base import array_of

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
