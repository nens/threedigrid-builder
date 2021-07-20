from threedigrid_builder.base import array_of
from threedigrid_builder.constants import CrossSectionShape

import numpy as np
import pygeos


__all__ = ["CrossSectionDefinitions"]


class CrossSectionDefinition:
    id: int
    code: str
    shape: CrossSectionShape
    height: float
    width: float


@array_of(CrossSectionDefinition)
class CrossSectionDefinitions:
    def compute_table(self):
        pass
