import numpy as np
cimport numpy as np

def set_2d_computational_nodes(
    double[:] origin,
    int kmax,
    int[:] mmax,
    int[:] nmax,
    double[:] dx,
    int[:,:] lg,
    int[:] nodk,
    int[:] nodm,
    int[:] nodn,
    int[:,:] quad_nod,
    double[:,:] bounds,
    double[:,:] coords,
):

    cdef int size_i = lg.shape[0]
    cdef int size_j = lg.shape[1]
    cdef int size_n = nodk.shape[0]

    f_set_2d_computational_nodes(
        origin=&origin[0],
        kmax=&kmax,
        mmax=&mmax[0],
        nmax=&nmax[0],
        dx=&dx[0],
        lg=&lg[0,0],
        size_i=&size_i,
        size_j=&size_j,
        nodk=&nodk[0],
        nodm=&nodm[0],
        nodn=&nodn[0],
        quad_nod=&quad_nod[0,0],
        bounds=&bounds[0,0],
        coords=&coords[0,0],
        size_n=&size_n
    )