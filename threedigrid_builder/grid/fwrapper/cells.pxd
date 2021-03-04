cdef extern:

    void* f_set_2d_computational_nodes(
        double *origin,
        int *kmax,
        int *mmax,
        int *nmax,
        double *dx,
        int *size_i,
        int *size_j,
        int *nodk,
        int *nodm,
        int *nodn,
        int *quad_idx,
        double *bounds,
        double *coords,
        int *size_n
    )

    void* f_set_2d_computational_lines(
        int *nodk,
        int *nodm,
        int *nodn,
        int *size_m,
        int *lgrmin,
        int *model_area,
        int *size_i,
        int *size_j,
        int *quad_idx,
        int *size_k,
        int *size_l,
        int *line,
        int *size_n,
    )