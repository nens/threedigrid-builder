cdef extern:

    void* init_nodes(int *handle, int *quadtree_handle)

    void* get_node_bnds(int *handle, int *n0, int *n1)
