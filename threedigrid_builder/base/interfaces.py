from .lines import Lines
from .nodes import Nodes
from abc import ABC
from abc import abstractmethod


__all__ = ["OutputInterface"]


class OutputInterface(ABC):
    """The metaclass (class template) for data output"""

    def __init__(self, path: str):
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