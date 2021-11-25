from .lines import Lines
from .nodes import Nodes
from .pumps import Pumps
from abc import ABC
from abc import abstractmethod
from pathlib import Path
from typing import Optional


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

    def __init__(self, path: Path, model_area_path: Optional[Path] = None):
        self.path = path
        self.model_area_path = model_area_path
        super().__init__()

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
