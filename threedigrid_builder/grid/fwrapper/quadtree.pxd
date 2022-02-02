cdef extern:

    void* make_quadtree(
        int *kmax,
        int *mmax,
        int *nmax,
        int *lgrmin,
        int *use_2d_flow,
        short *area_mask,
        int *lg,
        int *quad_idx,
        int *n0,
        int *n1,
        int *i0,
        int *i1,
        int *n_cells,
        int *n_line_u,
        int *n_line_v
    )