cdef extern:

    void* f_set_refinement(
        int *refine_id,
        double *refine_geom,
        int *n0,
        int *n1,
        int *refine_level,
        int *refine_type,
        double *bbox,
        int *mmax,
        int *nmax,
        double *dx,
        int *j0,
        int *lg,
        int *i0,
        int *i1
    )

    void* make_quadtree(
        int *kmax,
        int *mmax,
        int *nmax,
        int *lgrmin,
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