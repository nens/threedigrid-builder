cdef extern:

    void* f_set_2d_computational_nodes(
        double *origin,
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
        int *quad_nod,
        double *bounds,
        double *coords,
        int *size_n
    )