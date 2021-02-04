cimport quadtree
import numpy as np
cimport numpy as np
from enum import Enum
from enum import unique

@unique
class RefineType(Enum):
    LINESTRING = 1
    POLYGON = 2

cdef class QuadTree:
    cdef public int[2] handle

    def __init__(self, subgrid, int max_ref_lvl, min_grid_size):
        cdef double x0p = subgrid.origin[0]
        cdef double y0p = subgrid.origin[1]
        cdef double dxp = subgrid.pixel_size
        cdef int imax = subgrid.width
        cdef int jmax = subgrid.height
        cdef int min_pix_cell

        min_pix_cell = min_grid_size / subgrid.pixel_size

        quadtree.init_quadtree(
            handle=&self.handle[0],
            x0p=&x0p,
            y0p=&y0p,
            dxp=&dxp,
            imax=&imax,
            jmax=&jmax,
            lgrmin=&min_pix_cell,
            kmax=&max_ref_lvl
        )
        print("Handle: ", self.handle)

    cdef set_refinement(self, int id, double[::1,:] geom, int level, int type):
        cdef int status
        cdef int[:] c_handle = self.handle
        cdef int n0 = geom.shape[0]
        cdef int n1 = geom.shape[1]

        status = quadtree.set_refinement(
            handle=&c_handle[0],
            refine_id=&id,
            refine_geom=&geom[0, 0],
            n0=&n0,
            n1=&n1,
            refine_level=&level,
            refine_type=&type
        )
        return status            
    
    def set_refinements(self, refinements):
        cdef int refine_id 
        cdef double[::1,:] refine_geom
        cdef int[:] refine_level = refinements['refinement_level']
        cdef int refine_type
        cdef int i
        for i in range(len(refinements)):
            refine_id = refinements[i]['id']
            refine_geom = refinements[i]['geometry']
            refine_level = 
            refine_type = RefineType[refinements[i]['geometry_type'].upper()].value
            self.set_refinement(
                id = refine_id,
                geom = refine_geom,
                level = refine_level,
                type = refine_type
            )

    def create_active_grid(self, model_area_pix):
        cdef int[:] c_handle = self.handle
        cdef int[::1,:] model_area = model_area_pix
        cdef int n0 = model_area_pix.shape[0]
        cdef int n1 = model_area_pix.shape[1]
        quadtree.make_quadtree(handle=&c_handle[0])
        quadtree.set_active_2d_comp_cells(
            handle=&c_handle[0],
            model_area=&model_area[0,0],
            n0=&n0,
            n1=&n1
        )

    #def __del__(self):
