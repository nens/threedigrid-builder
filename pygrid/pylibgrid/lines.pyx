cimport lines
import numpy as np
cimport numpy as np


cdef class Lines:
    cdef public int[2] handle

    def __init__(self, nodes, model_area_pix):
        cdef int[:] c_nodes_handle = np.array(nodes.handle, dtype=np.int32)
        cdef int[::1,:] model_area = model_area_pix
        cdef int n0 = model_area_pix.shape[0]
        cdef int n1 = model_area_pix.shape[1]

        lines.init_lines(
            handle=&self.handle[0],
            nodes_handle=&c_nodes_handle[0],
            model_area=&model_area[0,0],
            n0=&n0,
            n1=&n1
        )