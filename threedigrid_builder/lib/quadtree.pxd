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

    int* make_quadtree(
        int *kmax,
        int *mmax,
        int *nmax,
        int *lgrmin,
        int *model_area,
        int *lg,
        int *m0,
        int *n0,
        int *n1,
        int *i0,
        int *i1
    )