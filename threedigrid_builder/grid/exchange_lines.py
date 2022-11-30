import pygeos

from threedigrid_builder.base import Array

__all__ = ["ExchangeLines"]


class ExchangeLine:
    id: int
    the_geom: pygeos.Geometry
    channel_id: int


class ExchangeLines(Array[ExchangeLine]):
    pass
