from . cimport quadtree
import numpy as np
cimport numpy as np
from enum import Enum
from enum import unique
import pygeos

@unique
class RefineType(Enum):
    LINESTRING = 1
    POLYGON = 3

cpdef set_refinement(
    int id,
    double[::1,:] geom,
    int level,
    int type,
    double[:] origin,
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

    status = quadtree.f_set_refinement(
        refine_id=&id,
        refine_geom=&geom[0, 0],
        n0=&n0,
        n1=&n1,
        refine_level=&level,
        refine_type=&type,
        origin=&origin[0],
        bbox=&bbox[0],
        mmax=&mmax[0],
        nmax=&nmax[0],
        dx=&dx[0],
        j0=&j0,
        lg=&lg[0,0],
        i0=&i0,
        i1=&i1
    )
            
    
def set_refinements(self, refinements, mmax, nmax, dx, lg):
    cdef int refine_id 
    cdef double[::1,:] refine_geom
    cdef int refine_level
    cdef int refine_type
    cdef int i
    pass


def create_quadtree(int kmax, int[:] mmax, int[:] nmax, int[:,:] lg):
    cdef int n0 = mmax.shape[0]
    cdef int i0 = lg.shape[0]
    cdef int i1 = lg.shape[1]
    quadtree.make_quadtree(
        kmax=&kmax,
        mmax=&mmax[0],
        nmax=&nmax[0],
        lg=&lg[0,0],
        n0=&n0,
        i0=&i0,
        i1=&i1
    )
    # quadtree.set_active_2d_comp_cells(
        # handle=&c_handle[0],
        # model_area=&model_area[0,0],
        # n0=&n0,
        # n1=&n1
    # )
