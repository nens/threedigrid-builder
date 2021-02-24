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
        int *lg,
        int *n0,
        int *i0,
        int *i1
    )

    #void set_active_2d_comp_cells(int *handle, int *model_area, int *n0, int *n1)
