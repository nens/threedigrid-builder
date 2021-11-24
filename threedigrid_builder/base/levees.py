from threedigrid_builder.base import array_of
from threedigrid_builder.constants import ContentType
from threedigrid_builder.constants import Material
from typing import Tuple

import numpy as np
import pygeos


__all__ = ["Levees", "Breaches"]


class Breach:
    id: int
    levl: int  # refers to line
    levbr: float
    levmat: Material
    content_pk: int  # refers to v2_connected_pnt
    coordinates: Tuple[float, float]


@array_of(Breach)
class Breaches:
    pass


class Levee:
    id: int
    the_geom: pygeos.Geometry
    crest_level: float
    max_breach_depth: float
    material: Material


@array_of(Levee)
class Levees:
    def get_breaches(self, nodes, lines, connected_points):
        # get lines that have a connected point
        (line_idx,) = np.where(
            lines.content_type == ContentType.TYPE_V2_ADDED_CALCULATION_POINT
        )
        # refine to those connected points that have a levee
        conn_pnt_idx = connected_points.id_to_index(lines.content_pk[line_idx])
        has_levee = connected_points.levee_id[conn_pnt_idx] != -9999

        line_idx = line_idx[has_levee]
        conn_pnt_idx = conn_pnt_idx[has_levee]
        levee_idx = self.id_to_index(connected_points.levee_id[conn_pnt_idx])

        if not np.any(has_levee):
            return Breaches(id=[])

        # set levee crest level (if connected point crest level is undefined)
        mask = np.isnan(connected_points.exchange_level[conn_pnt_idx]) & np.isfinite(
            self.crest_level[levee_idx]
        )
        lines.dpumax[line_idx[mask]] = self.crest_level[levee_idx[mask]]

        # only consider levees that have the correct properties set
        mask = (np.isfinite(self.max_breach_depth) & (self.material != -9999))[
            levee_idx
        ]
        if not np.any(mask):
            return Breaches(id=[])

        line_idx = line_idx[mask]
        levee_idx = levee_idx[mask]

        # compute the intersections
        lines_breaches = lines[line_idx]
        lines_breaches.set_line_coords(nodes)
        lines_breaches.line_geometries[:] = None
        lines_breaches.fix_line_geometries()
        points = pygeos.intersection(
            lines_breaches.line_geometries, self.the_geom[levee_idx]
        )
        # edge case: non-point type intersections
        non_point = pygeos.get_type_id(points) != 0
        points[non_point] = pygeos.point_on_surface(points[non_point])

        return Breaches(
            id=range(len(line_idx)),
            levl=lines.index_to_id(line_idx),
            levbr=self.max_breach_depth[levee_idx],
            levmat=self.material[levee_idx],
            content_pk=lines.content_pk[line_idx],
            coordinates=np.array([pygeos.get_x(points), pygeos.get_y(points)]).T,
        )
