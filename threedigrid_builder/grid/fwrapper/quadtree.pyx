import numpy as np

cimport numpy as np

from enum import Enum
from enum import unique


@unique
class RefineType(Enum):
    LINESTRING = 1
    POLYGON = 3

cpdef set_refinement(
    int id,
    double[::1,:] geom,
    int level,
    int type,
    double[:] bbox,
    int[:] mmax,
    int[:] nmax,
    double[:] dx,
    int[::1,:] lg
    ):

    cdef int n0 = geom.shape[0]
    cdef int n1 = geom.shape[1]
    cdef int j0 = mmax.shape[0]
    cdef int i0 = lg.shape[0]
    cdef int i1 = lg.shape[1]

    status = f_set_refinement(
        refine_id=&id,
        refine_geom=&geom[0, 0],
        n0=&n0,
        n1=&n1,
        refine_level=&level,
        refine_type=&type,
        bbox=&bbox[0],
        mmax=&mmax[0],
        nmax=&nmax[0],
        dx=&dx[0],
        j0=&j0,
        lg=&lg[0,0],
        i0=&i0,
        i1=&i1
    )

cpdef create_quadtree(
    int kmax, 
    int[:] mmax,
    int[:] nmax,
    int lgrmin,
    short[::1,:] area_mask,
    int[::1,:] quad_idx,
    int[::1,:] lg,
):

    cdef int n0 = area_mask.shape[0]
    cdef int n1 = area_mask.shape[1]
    cdef int i0 = lg.shape[0]
    cdef int i1 = lg.shape[1]
    cdef int n_cells
    cdef int n_line_u
    cdef int n_line_v

    make_quadtree(
        kmax=&kmax,
        mmax=&mmax[0],
        nmax=&nmax[0],
        lgrmin=&lgrmin,
        area_mask=&area_mask[0,0],
        lg=&lg[0,0],
        quad_idx=&quad_idx[0,0],
        n0=&n0,
        n1=&n1,
        i0=&i0,
        i1=&i1,
        n_cells=&n_cells,
        n_line_u=&n_line_u,
        n_line_v=&n_line_v
    )

    return <object>n_cells, <object>(n_line_u, n_line_v)
