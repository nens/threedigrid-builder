from threedigrid_builder.base import Lines
from threedigrid_builder.base import Nodes
from .quadtree import QuadTree


__all__ = ["Grid"]


class Grid:
    def __init__(self, nodes: Nodes, lines: Lines):
        if not isinstance(nodes, Nodes):
            raise TypeError(f"Expected Nodes instance, got {type(nodes)}")
        if not isinstance(lines, Lines):
            raise TypeError(f"Expected Lines instance, got {type(lines)}")
        self.nodes = nodes
        self.lines = lines

    def concatenate(self, other):
        """Concatenate two grids without renumbering nodes."""
        if self.__class__ is not other.__class__:
            raise TypeError(
                "Cannot concatenate {self} with {other} as they are not of "
                "equal types."
            )
        return self.__class__(
            nodes=self.nodes + other.nodes,
            lines=self.lines + other.lines,
        )

    def __repr__(self):
        return f"<Grid object with {len(self.nodes)} nodes and {len(self.lines)} lines>"

    @classmethod
    def from_quadtree(
        cls, subgrid_meta, num_refine_levels, min_gridsize, refinements
    ):
        """Construct the 2D grid based on the quadtree object.
        """
        quadtree = QuadTree(
            subgrid_meta,
            num_refine_levels,
            min_gridsize,
            refinements
        )

        nodes = quadtree.get_nodes(subgrid_meta)
        lines = Lines(id=[])

        return cls(nodes=nodes, lines=lines)
