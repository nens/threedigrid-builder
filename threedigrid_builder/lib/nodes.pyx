cimport nodes
cimport quadtree
import numpy as np
from numpy cimport npy_intp, ndarray, NPY_DOUBLE, NPY_INT
from utils cimport create_array


cdef class Nodes:
    cdef public int[2] handle

    def __init__(self, quadtree):
        cdef int[:] c_quadtree_handle = np.array(quadtree.handle, dtype=np.int32)
        nodes.init_nodes(
            handle=&self.handle[0],
            quadtree_handle=&c_quadtree_handle[0]
        )

    @property
    def nodtot(self):
        return <object>nodes.get_nodtot(handle=&self.handle[0])

    def get_coordinates(self, int n0, int n1):
        cdef rank = 2
        cdef void* data_ptr
        cdef npy_intp shape[2]
        cdef ndarray arr
        shape[0] = n1 - n0 + 1
        shape[1] = 2

        data_ptr = nodes.get_coordinates(
            handle=&self.handle[0],
            n0=&n0,
            n1=&n1
        )
        arr = create_array(rank, shape, NPY_DOUBLE, data_ptr)
        return <object>arr

    def get_bnds(self, int n0, int n1):
        cdef rank = 2
        cdef void* data_ptr
        cdef npy_intp shape[2]
        cdef ndarray arr
        shape[0] = n1 - n0 + 1
        shape[1] = 4

        data_ptr = nodes.get_node_bnds(
            handle=&self.handle[0],
            n0=&n0,
            n1=&n1
        )
        arr = create_array(rank, shape, NPY_DOUBLE, data_ptr)
        return <object>arr

    def get_id(self, int n0, int n1):
        cdef rank = 1
        cdef void* data_ptr
        cdef npy_intp shape[1]
        shape[0] = n1 - n0 + 1

        data_ptr = nodes.get_node_id(
            handle=&self.handle[0],
            n0=&n0,
            n1=&n1
        )
        arr = create_array(rank, shape, NPY_INT, data_ptr)
        return <object>arr
