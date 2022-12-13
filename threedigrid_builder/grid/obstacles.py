import numpy as np
import pygeos

from threedigrid_builder.base import Array, Lines


class Obstacle:
    id: int
    crest_level: float
    the_geom: pygeos.Geometry


class Obstacles(Array[Obstacle]):
    def compute_dpumax(self, lines: Lines, where):
        """Compute the dpumax for lines that intersect with the obstacles.

        Args:
            lines (Lines)
            where (int array): indices into lines to compute dpumax for

        Returns:
            dpumax for each line (length equals length of 'where' array)
        """
        exchange_level = np.full(len(where), np.nan)
        if len(self) == 0:
            return exchange_level
        coordinates = lines.line_coords[where]
        if not np.isfinite(coordinates).all():
            raise ValueError(
                f"{lines.__class__.__name__} object has inclomplete line_coords."
            )
        lines_tree = pygeos.STRtree(pygeos.linestrings(coordinates.reshape(-1, 2, 2)))
        inscts = lines_tree.query_bulk(self.the_geom, predicate="intersects")
        for i in range(len(self)):
            indices = inscts[1, np.where(inscts[0, :] == i)]
            exchange_level[indices] = np.fmax(
                exchange_level[indices], self.crest_level[i]
            )
        return exchange_level
