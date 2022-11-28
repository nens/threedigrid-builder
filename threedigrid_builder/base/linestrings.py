import pygeos
from .array import array_of
import numpy as np

__all__ = ["LineStrings", "PointsOnLine"]

class LineString:
    id: int
    the_geom: pygeos.Geometry  # LineString

    
class PointOnLine:
    """A point determined by its position 'along' a linestring"""
    id: int  # just a unique number, mostly unused
    content_pk: int  # externally determined id (e.g. cs location id)
    s1d: float  # the position along the linestring
    s1d_cum: float
    linestring_idx: int


@array_of(LineString)
class LineStrings:
    @property
    def length(self):
        return pygeos.length(self.the_geom)
    
    @property
    def cum_length(self):
        result = np.zeros(len(self) + 1, dtype=float)
        np.cumsum(self.length + 1e-6, out=result[1:])
        return result

    def locate_points(self, geoms, ids, linestring_ids) -> "PointsOnLine":
        return PointsOnLine.from_geometries(self, geoms, ids, linestring_ids)

    
@array_of(PointOnLine)
class PointsOnLine:
    @classmethod
    def from_geometries(cls, linestrings: LineStrings, geoms, linestring_ids, **kwargs):
        linestring_idx = linestrings.id_to_index(linestring_ids, check_exists=True)
        s1d = pygeos.line_locate_point(linestrings.the_geom[linestring_idx], geoms)
        return cls.from_s1d(linestrings, s1d, linestring_ids, **kwargs)
    
    @classmethod
    def from_s1d(cls, linestrings: LineStrings, s1d, linestring_ids, **kwargs):
        if np.any(~np.isfinite(s1d)):
            raise ValueError("NaN values encountered in s1d")
        linestring_idx = linestrings.id_to_index(linestring_ids, check_exists=True)
        result = cls(
            id=np.arange(len(s1d)),
            s1d=s1d,
            s1d_cum=s1d + linestrings.cum_length[linestring_idx],
            linestring_idx=linestring_idx,
            **kwargs,
        )
        result.reorder_by("s1d_cum")
        return result
    
    def neighbors(self, other: "PointOnLine"):
        """Return the (indices of) the points that are before and after self.
                
        'other' must contain at least 1 point per linestring.

        If a point does not have a neigbour before it, the two returned indices
        will be equal (both pointing to the neighbour after it). Vice versa, if it
        does not have a neighboar after it, the two returned indices will be that
        of the neighboar before.
        """
        idx_2 = np.searchsorted(other.s1d_cum, self.s1d_cum)
        idx_1 = idx_2 - 1
        out_of_bounds_1 = idx_1 < 0
        out_of_bounds_2 = idx_2 >= len(other)
        idx_1[out_of_bounds_1] = 0
        idx_2[out_of_bounds_2] = len(other) - 1

        # Fix situations where idx_1 is incorrect
        extrap_mask = (other.linestring_idx[idx_1] != self.linestring_idx) | out_of_bounds_1
        idx_1[extrap_mask] = idx_2[extrap_mask]
        idx_2[extrap_mask] = np.clip(idx_2[extrap_mask] + 1, None, len(other) - 1)
        equalize = extrap_mask & (other.linestring_idx[idx_2] != self.linestring_idx)
        idx_2[equalize] -= 1

        # Fix situations where idx_2 is incorrect
        extrap_mask = (other.linestring_idx[idx_2] != self.linestring_idx) | out_of_bounds_2
        idx_2[extrap_mask] = idx_1[extrap_mask]
        idx_1[extrap_mask] = np.clip(idx_2[extrap_mask] - 1, 0, None)
        equalize = extrap_mask & (other.linestring_idx[idx_1] != self.linestring_idx)
        idx_1[equalize] += 1

        return idx_1, idx_2
