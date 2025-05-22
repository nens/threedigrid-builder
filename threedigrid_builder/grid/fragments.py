import shapely

from threedigrid_builder.base import Array

__all__ = ["Fragment", "Fragments"]


class Fragment:
    id: int
    node_id: int
    the_geom: shapely.Geometry


class Fragments(Array[Fragment]):
    pass
