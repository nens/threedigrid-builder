cimport lines
import numpy as np
from numpy cimport npy_intp, ndarray, NPY_DOUBLE, NPY_INT
from utils cimport create_array


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
    
    @property
    def lintot(self):
        return <object>lines.get_lintot(handle=&self.handle[0])

    def get_line(self, int n0, int n1):
        cdef rank = 2
        cdef void* data_ptr
        cdef npy_intp shape[2]
        cdef ndarray arr
        shape[0] = n1 - n0 + 1
        shape[1] = 2

        data_ptr = lines.get_line(
            handle=&self.handle[0],
            n0=&n0,
            n1=&n1
        )
        arr = create_array(rank, shape, NPY_INT, data_ptr)
        return <object>arr

    def get_id(self, int n0, int n1):
        cdef rank = 1
        cdef void* data_ptr
        cdef npy_intp shape[1]
        shape[0] = n1 - n0 + 1

        data_ptr = lines.get_line_id(
            handle=&self.handle[0],
            n0=&n0,
            n1=&n1
        )
        arr = create_array(rank, shape, NPY_INT, data_ptr)
        return <object>arr