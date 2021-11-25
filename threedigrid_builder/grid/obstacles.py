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


def apply_obstacles(lines, obstacles, levees):
    """Set obstacles on 2D lines by calculating intersection between
    2D flowline line_coords and obstacle linestring. In case of an
    intersection the flowline becomes and LINE_2D_OBSTACLE.
    The LINE_2D_OBSTACLE is set on kcu. Also flod and flou of lines is set to
    crest_level.

    Args:
        lines (Lines)
        obstacles (Obstacles)
        levees (Levees)
    """
    if len(obstacles) == 0 and len(levees) == 0:
        return
    the_geom = np.concatenate([obstacles.the_geom, levees.the_geom])
    crest_level = np.concatenate([obstacles.crest_level, levees.crest_level])

    is_2d = np.isin(lines.kcu, (LineType.LINE_2D_U, LineType.LINE_2D_V))
    coordinates = lines.line_coords[is_2d]
    if not np.isfinite(coordinates).all():
        raise ValueError(
            f"{lines.__class__.__name__} object has inclomplete line_coords."
        )
    lines_tree = pygeos.STRtree(pygeos.linestrings(coordinates.reshape(-1, 2, 2)))

    inscts = lines_tree.query_bulk(the_geom, predicate="intersects")
    is_u = np.where(lines.kcu == LineType.LINE_2D_U)[0]
    mask = np.isin(is_u, inscts[1, :])
    lines.kcu[is_u[mask]] = LineType.LINE_2D_OBSTACLE_U
    is_v = np.where(lines.kcu == LineType.LINE_2D_V)[0]
    mask = np.isin(is_v, inscts[1, :])
    lines.kcu[is_v[mask]] = LineType.LINE_2D_OBSTACLE_V
    for i in range(len(crest_level)):
        indices = inscts[1, np.where(inscts[0, :] == i)]
        lines.flod[indices] = np.fmax(lines.flod[indices], crest_level[i])
        lines.flou[indices] = lines.flod[indices]
