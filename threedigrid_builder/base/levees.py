from threedigrid_builder.base import array_of
from threedigrid_builder.constants import ContentType, LineType
from threedigrid_builder.constants import Material
from typing import Tuple

import numpy as np
import pygeos


__all__ = ["Levees", "Breaches"]


class Breach:
    id: int
    line_id: int
    content_pk: int  # refers to v2_connected_pnt
    coordinates: Tuple[float, float]
    levee_id: int


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
        points[pygeos.is_empty(points)] = None

        return Breaches(
            id=range(len(line_idx)),
            line_id=lines.index_to_id(line_idx),
            levee_id=self.index_to_id(levee_idx),
            content_pk=lines.content_pk[line_idx],
            coordinates=np.array([pygeos.get_x(points), pygeos.get_y(points)]).T,
        )

    def get_breaches2(self, nodes, lines):
        # only consider levees that have the correct properties set
        mask = (np.isfinite(self.max_breach_depth) & (self.material != -9999))
        if not np.any(mask):
            return Breaches(id=[])

        # get 1D2D lines that may cross a levee
        LINETYPE_1D2D_OPEN = [
            LineType.LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER,
            LineType.LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER
        ]
        lines_1d2d_open = lines[np.isin(lines.kcu, LINETYPE_1D2D_OPEN)]
        if len(lines_1d2d_open) == 0:
            return Breaches(id=[])

        # get crossings
        lines_1d2d_open.set_line_coords(nodes)
        lines_1d2d_open.fix_line_geometries()
        tree = pygeos.STRtree(lines_1d2d_open.line_geometries)

        levee_idx, line_idx = tree.query_bulk(self.the_geom[mask], predicate="crosses")
        levee_idx = np.where(mask)[0][levee_idx]

        # take the first levee for each line
        line_idx, where = np.unique(line_idx, return_index=True)
        levee_idx = levee_idx[where]

        # compute the intersections
        points = pygeos.intersection(
            lines_1d2d_open.line_geometries[line_idx], self.the_geom[levee_idx]
        )
        # edge case: non-point type intersections
        non_point = pygeos.get_type_id(points) != 0
        points[non_point] = pygeos.point_on_surface(points[non_point])
        points[pygeos.is_empty(points)] = None

        return Breaches(
            id=range(len(line_idx)),
            line_id=lines_1d2d_open.index_to_id(line_idx),
            levee_id=self.index_to_id(levee_idx),
            coordinates=np.array([pygeos.get_x(points), pygeos.get_y(points)]).T,
        )
