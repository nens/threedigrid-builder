import numpy as np
import shapely

from threedigrid_builder.base import Array, Lines


class Obstacle:
    id: int
    crest_level: float
    the_geom: shapely.Geometry


class Obstacles(Array[Obstacle]):
    def compute_dpumax(self, lines: Lines, where):
        """Compute the dpumax for lines that intersect with the obstacles.

        Args:
            lines (Lines)
            where (int array): indices into lines to compute dpumax for

        Returns:
            dpumax for each line (length equals length of 'where' array)
        """
        obstacle_idx = np.full(len(where), -9999, dtype=np.int32)
        exchange_level = np.full(len(where), np.nan)
        if len(self) == 0:
            return exchange_level, obstacle_idx
        line_geometries = lines.line_geometries[where]
        if shapely.is_missing(line_geometries).any():
            raise ValueError(
                "Obstacles.compute_dpumax requires line_geometries to be set"
            )
        lines_tree = shapely.STRtree(lines.line_geometries[where])
        inscts = lines_tree.query(self.the_geom, predicate="intersects")
        for i in range(len(self)):
            line_idx = inscts[1, np.where(inscts[0, :] == i)[0]]
            is_greater = ~(self.crest_level[i] <= exchange_level[line_idx])
            exchange_level[line_idx[is_greater]] = self.crest_level[i]
            obstacle_idx[line_idx[is_greater]] = i
        return exchange_level, obstacle_idx
