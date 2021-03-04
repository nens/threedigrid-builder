import numpy as np
cimport numpy as np

def set_2d_computational_nodes(
    double[:] origin,
    int kmax,
    int[:] mmax,
    int[:] nmax,
    double[:] dx,
    int[:] nodk,
    int[:] nodm,
    int[:] nodn,
    int[::1,:] quad_idx,
    double[::1,:] bounds,
    double[::1,:] coords,
):

    cdef int size_i = quad_idx.shape[0]
    cdef int size_j = quad_idx.shape[1]
    cdef int size_n = nodk.shape[0]

    f_set_2d_computational_nodes(
        origin=&origin[0],
        kmax=&kmax,
        mmax=&mmax[0],
        nmax=&nmax[0],
        dx=&dx[0],
        size_i=&size_i,
        size_j=&size_j,
        nodk=&nodk[0],
        nodm=&nodm[0],
        nodn=&nodn[0],
        quad_idx=&quad_idx[0,0],
        bounds=&bounds[0,0],
        coords=&coords[0,0],
        size_n=&size_n
    )


def set_2d_computational_lines(
    int[:] nodk,
    int[:] nodm,
    int[:] nodn,
    int lgrmin,
    int[::1,:] model_area,
    int[::1,:] quad_idx,
    int[::1,:] line,
):

    cdef int size_m = nodk.shape[0]
    cdef int size_i = model_area.shape[0]
    cdef int size_j = model_area.shape[1]
    cdef int size_k = quad_idx.shape[0]
    cdef int size_l = quad_idx.shape[1]
    cdef int size_n = line.shape[0]

    f_set_2d_computational_lines(
        nodk=&nodk[0],
        nodm=&nodm[0],
        nodn=&nodn[0],
        size_m=&size_m,
        lgrmin=&lgrmin,
        model_area=&model_area[0,0],
        size_i=&size_i,
        size_j=&size_j,
        quad_idx=&quad_idx[0,0],
        size_k=&size_k,
        size_l=&size_l,
        line=&line[0,0],
        size_n=&size_n
    )