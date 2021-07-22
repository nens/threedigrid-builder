from threedigrid_builder.base import array_of
from threedigrid_builder.constants import LineType

import numpy as np
import pygeos


class Obstacle:
    id: int
    crest_level: float
    the_geom: pygeos.Geometry


@array_of(Obstacle)
class Obstacles:
    pass


def apply_obstacles(lines, obstacles):
    """Set obstacles on 2D lines by calculating intersection between
    2D flowline line_coords and obstacle linestring. In case of an
    intersection the flowline becomes and LINE_2D_OBSTACLE.
    The LINE_2D_OBSTACLE is set on kcu. Also flod and flou of lines is set to
    crest_level.

    Args:
        lines (Lines)
        obstacles (Obstacles)
    """
    is_2d = np.isin(lines.kcu, (LineType.LINE_2D_U, LineType.LINE_2D_V))
    coordinates = lines.line_coords[is_2d]
    if not np.isfinite(coordinates).all():
        raise ValueError(f"{lines.__class__.__name__} object has inclomplete line_coords.")
    lines_tree = pygeos.STRtree(pygeos.linestrings(coordinates.reshape(-1, 2, 2)))

    inscts = lines_tree.query_bulk(obstacles.the_geom, predicate="intersects")
    lines.kcu[inscts[1, :]] = LineType.LINE_2D_OBSTACLE
    for i in range(len(obstacles.id)):
        indices = inscts[1, np.where(inscts[0, :] == i)]
        lines.flod[indices] = np.fmax(
            lines.flod[indices], obstacles.crest_level[i]
        )
        lines.flou[indices] = lines.flod[indices]
