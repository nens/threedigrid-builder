from threedigrid_builder.base import array_of
from threedigrid_builder.constants import PumpType
from typing import Tuple


__all__ = ["Pumps"]


class Pump:
    id: int
    code: str
    capacity: float
    connection_node_start_id: int
    connection_node_end_id: int  # optional!
    type: PumpType
    start_level: float
    lower_stop_level: float
    upper_stop_level: float
    # zoom_category
    # display_name
    # sewerage
    bottom_level: float
    content_pk: int
    line: Tuple[int, int]
    line_coords: Tuple[float, float, float, float]


@array_of(Pump)
class Pumps:
    pass
