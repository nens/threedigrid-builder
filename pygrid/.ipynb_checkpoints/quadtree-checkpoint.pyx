cimport quadtree
import numpy as np

class Quadtree:

    def __init__(self, double x0p, double y0p, double dxp, int imax, int jmax, int lgrmin, int kmax):
        # cdef double *c_x0p
        # cdef double *c_y0p
        # cdef double *c_dxp
        # cdef int *c_imax
        # cdef int *c_jmax
        # cdef int *c_lgrmin
        # cdef int *c_kmax
        cdef int result[2]
        cdef int[:] handle = np.zeros((2,), dtype=np.int32)
        quadtree.init_quadtree(handle=result, x0p=&x0p, y0p=&y0p, dxp=&dxp, imax=&imax, jmax=&jmax, lgrmin=&lgrmin, kmax=&kmax)
        handle = result
        #self._handle = handle
        print('RESULT: ', handle)
        #self._handle = handle[0] if isinstance(handle, tuple) else handle
        return
    