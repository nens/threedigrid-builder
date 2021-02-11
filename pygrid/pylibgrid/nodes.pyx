from numpy cimport ndarray, import_array, npy_intp, PyArray_SimpleNewFromData, \
                   NPY_DOUBLE, NPY_F_CONTIGUOUS
cimport nodes
cimport quadtree
import numpy as np
from cpython cimport PyObject
from cpython.ref cimport PyTypeObject

# Initialization
import_array()

cdef class Nodes:
    cdef public int[2] handle

    def __init__(self, quadtree):
        cdef int[:] c_quadtree_handle = np.array(quadtree.handle, dtype=np.int32)
        nodes.init_nodes(
            handle=&self.handle[0],
            quadtree_handle=&c_quadtree_handle[0]
        )

    def node_bnds(self, int n0, int n1):
        cdef void* data_ptr
        cdef npy_intp shape[2]
        cdef ndarray ndarray
        shape[1] = n1 - n0 + 1
        shape[0] = 4

        data_ptr = nodes.get_node_bnds(
            handle=&self.handle[0],
            n0=&n0,
            n1=&n1
        )

        ndarray = PyArray_SimpleNewFromData(2, shape, NPY_DOUBLE, data_ptr)

        # ndarray = PyArray_NewFromDescr(
            # PyArray_Type,
            # PyArray_DescrFromType(NPY_DOUBLE),
            # 2,
            # shape,
            # NULL,
            # data_ptr,
            # NPY_F_CONTIGUOUS,
            # None
        # )

        return <object>ndarray