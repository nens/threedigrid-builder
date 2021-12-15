cdef extern:

    void* f_set_2d_computational_nodes_lines(
        double *origin,
        int *lgrmin,
        int *kmax,
        int *mmax,
        int *nmax,
        double *dx,
        int *lg,
        int *size_i,
        int *size_j,
        int *nodk,
        int *nodm,
        int *nodn,
        int *quad_idx,
        double *bounds,
        double *coords,
        int *pixel_coords,
        int *size_n,
        short *area_mask,
        int *size_a,
        int *size_b,
        int *line,
        int *cross_pix_coords,
        int *n_line_u,
        int *n_line_v
    )
