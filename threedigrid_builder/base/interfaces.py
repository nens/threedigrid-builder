from .lines import Lines
from .nodes import Nodes
from .pumps import Pumps
from abc import ABC
from abc import abstractmethod
from pathlib import Path


__all__ = ["OutputInterface"]


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
