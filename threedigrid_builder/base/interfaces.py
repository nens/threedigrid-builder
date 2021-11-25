from .lines import Lines
from .nodes import Nodes
from .pumps import Pumps
from abc import ABC
from abc import abstractmethod
from pathlib import Path
from typing import Optional
from threedigrid_builder.exceptions import SchematisationError


__all__ = ["OutputInterface", "RasterInterface"]


class OutputInterface(ABC):
    """The metaclass (class template) for data output"""

    def __init__(self, path: Path):
        self.path = path
        super().__init__()

    @abstractmethod
    def __enter__(self):
        pass

    @abstractmethod
    def __exit__(self, type, value, traceback):
        pass

    @abstractmethod
    def write_nodes(self, nodes: Nodes):
        pass

    @abstractmethod
    def write_lines(self, lines: Lines):
        pass

    @abstractmethod
    def write_pumps(self, pumps: Pumps):
        pass

    @abstractmethod
    def write_cross_sections(self, cross_sections):
        pass

    @abstractmethod
    def write_nodes_embedded(self, nodes_embedded: Nodes):
        pass


class RasterInterface(ABC):
    """The metaclass (class template) for raster data input"""
    GT_TOLERANCE = 7

    def __init__(self, path: Path, model_area_path: Optional[Path] = None, epsg_code: Optional[int] = None):
        self.path = path
        self.model_area_path = model_area_path
        self.transform = None
        self.epsg_code = None
        super().__init__()

    def set_transform(self, transform):
        """Transform should be given in RasterIO order:

        >>> dx, 0, x, 0, dy, y
        """
        a, b, c, d, e, f = [round(x, self.GT_TOLERANCE) for x in transform]

        # check if pixels are square and if cross terms are 0
        if abs(a) != abs(e):
            raise SchematisationError(f"Raster pixels are non-square ({abs(a)}x{abs(e)}).")
        if b != 0.0 or d != 0.0:
            raise SchematisationError(f"Raster pixel grid has tilt (cross terms: {b},{d}).")

        self.transform = (a, b, c, d, e, f)

    def set_epsg_code(self, epsg_code):
        if epsg_code is None:
            return
        if self.epsg_code is not None and epsg_code != self.epsg_code:
            raise SchematisationError(f"Raster and SQLite epsg code mismatch ({epsg_code} != {self.epsg_code}).")
        self.epsg_code = epsg_code

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
