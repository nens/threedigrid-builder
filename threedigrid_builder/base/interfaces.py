from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional

from pyproj import CRS

from threedigrid_builder.exceptions import SchematisationError

__all__ = ["OutputInterface", "RasterInterface"]


class OutputInterface(ABC):
    """The metaclass (class template) for data output"""

    def __init__(self, path: Path):
        self.path = path
        super().__init__()

    @staticmethod
    def available():
        return True

    @abstractmethod
    def __enter__(self):
        pass

    @abstractmethod
    def __exit__(self, type, value, traceback):
        pass

    @abstractmethod
    def write(self, grid):
        pass


class RasterInterface(ABC):
    """The metaclass (class template) for raster data input"""

    GT_TOLERANCE = 7

    def __init__(
        self,
        path: Path,
        model_area_path: Optional[Path] = None,
    ):
        self.path = path
        self.model_area_path = model_area_path
        self.transform = None
        self.crs = None
        super().__init__()

    def set_transform(self, transform):
        """Transform should be given in RasterIO order:

        >>> dx, 0, x, 0, dy, y
        """
        a, b, c, d, e, f = [round(x, self.GT_TOLERANCE) for x in transform]

        # check if pixels are square and if cross terms are 0
        if abs(a) != abs(e):
            raise SchematisationError(
                f"Raster pixels are non-square ({abs(a)}x{abs(e)})."
            )
        if b != 0.0 or d != 0.0:
            raise SchematisationError(
                f"Raster pixel grid has tilt (cross terms: {b},{d})."
            )

        self.transform = (a, b, c, d, e, f)

    def set_crs(self, wkt: str):
        self.crs = CRS.from_wkt(wkt)

    @property
    def pixel_size(self):
        return abs(self.transform[0])

    def __repr__(self):
        return f"<{self.__class__.__name__} of file {self.path}>"

    @abstractmethod
    def __enter__(self):
        pass

    @abstractmethod
    def __exit__(self, type, value, traceback):
        pass

    @abstractmethod
    def read(self):
        pass
