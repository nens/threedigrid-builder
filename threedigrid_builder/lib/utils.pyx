from numpy cimport ndarray, import_array, npy_intp, PyArray_SimpleNewFromData, \
                   PyArray_New, NPY_DOUBLE, NPY_INT, NPY_FLOAT64, NPY_F_CONTIGUOUS

from cpython cimport PyObject
from cpython.ref cimport PyTypeObject


# Initialization
import_array()


cdef create_array(int rank, npy_intp* shape, int dtype, void* data_ptr):
    arr = PyArray_New(ndarray, rank, shape, dtype, NULL, data_ptr, dtype, NPY_F_CONTIGUOUS, None)
    return arr