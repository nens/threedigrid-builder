from typing import Tuple

from .array import Array

__all__ = ["Quarters"]


class Quarter:
    id: int
    line: Tuple[int, int]
    neighbour_node: Tuple[int, int]


class Quarters(Array[Quarter]):
    """Calculation quarters."""
