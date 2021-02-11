from numpy cimport npy_intp

cdef object create_array(int rank, npy_intp* shape, int dtype, void* data_ptr)
