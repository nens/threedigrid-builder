import numpy as np

cimport numpy as np


def set_2d_computational_nodes_lines(
    double[:] origin,
    int lgrmin,
    int kmax,
    int[:] mmax,
    int[:] nmax,
    double[:] dx,
    int[::1,:] lg,
    int[::1,:] quad_idx,
    short[::1,:] area_mask,
    int[:] nodk,
    int[:] nodm,
    int[:] nodn,
    double[::1,:] bounds,
    double[::1,:] coords,
    int[::1,:] pixel_coords,
    int[::1,:] line,
    int[::1,:] cross_pix_coords,
    int n_line_u,
    int n_line_v
):

    cdef int size_i = quad_idx.shape[0]
    cdef int size_j = quad_idx.shape[1]
    cdef int size_n = nodk.shape[0]
    cdef int size_a = area_mask.shape[0]
    cdef int size_b = area_mask.shape[1]

    f_set_2d_computational_nodes_lines(
        origin=&origin[0],
        lgrmin=&lgrmin,
        kmax=&kmax,
        mmax=&mmax[0],
        nmax=&nmax[0],
        dx=&dx[0],
        lg=&lg[0, 0],
        size_i=&size_i,
        size_j=&size_j,
        nodk=&nodk[0],
        nodm=&nodm[0],
        nodn=&nodn[0],
        quad_idx=&quad_idx[0, 0],
        bounds=&bounds[0, 0],
        coords=&coords[0, 0],
        pixel_coords=&pixel_coords[0, 0],
        size_n=&size_n,
        area_mask=&area_mask[0, 0],
        size_a=&size_a,
        size_b=&size_b,
        line=&line[0, 0] if line.size > 0 else NULL,
        cross_pix_coords=&cross_pix_coords[0, 0] if cross_pix_coords.size > 0 else NULL,
        n_line_u=&n_line_u,
        n_line_v=&n_line_v
    )
