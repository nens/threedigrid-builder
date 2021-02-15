cdef extern:

    void init_quadtree(int *handle, double *x0p, double *y0p, double *dxp, int *imax, int *jmax, int *lgrmin, int *kmax)

    int set_refinement(int*handle, int *refine_id, double *refine_geom, int *n0, int *n1, int *refine_level, int *refine_type)

    void make_quadtree(int *handle)

    void set_active_2d_comp_cells(int *handle, int *model_area, int *n0, int *n1)
