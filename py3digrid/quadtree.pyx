cimport quadtree
import numpy as np
cimport numpy as np
from enum import Enum
from enum import unique

@unique
class RefineType(Enum):
    LINESTRING = 1
    POLYGON = 2

class QuadTree:

    def __init__(self, subgrid, int max_ref_lvl, min_grid_size):
        cdef double x0p = subgrid.origin[0]
        cdef double y0p = subgrid.origin[1]
        cdef double dxp = subgrid.pixel_size
        cdef int imax = subgrid.width
        cdef int jmax = subgrid.height,
        cdef int min_pix_cell,

        min_pix_cell = min_grid_size / subgrid.pixel_size

        self._handle = np.zeros(2, dtype=np.int32)
        cdef int[:] c_handle = self._handle
        quadtree.init_quadtree(handle=&c_handle[0], x0p=&x0p, y0p=&y0p, dxp=&dxp, imax=&imax, jmax=&jmax, lgrmin=&min_pix_cell, kmax=&max_ref_lvl)
        print("Handle: ", self._handle)

    def set_refinement(self, refinement):
        cdef int refine_id = refinement['id']
        cdef double[:,:] refine_geom = refinement['geometry'].T
        cdef int n0 = refine_geom.shape[0]
        cdef int n1 = refine_geom.shape[1]
        cdef int refine_level = refinement['refinement_level']
        cdef int refine_type = RefineType[refinement['geometry_type'].upper()].value
        cdef int status
        cdef int[:] c_handle = self._handle

        status = quadtree.set_refinement(handle=&c_handle[0],
                                refine_id=&refine_id,
                                refine_geom=&refine_geom[0, 0],
                                n0=&n0,
                                n1=&n1,
                                refine_level=&refine_level,
                                refine_type=&refine_type)
    
        return status
            

