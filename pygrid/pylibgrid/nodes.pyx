cimport nodes
cimport quadtree
import numpy as np
cimport numpy as np


cdef class Nodes:
    cdef public int[2] handle

    def __init__(self, quadtree):
        cdef int[:] c_quadtree_handle = np.array(quadtree.handle, dtype=np.int32)
        nodes.init_nodes(
            handle=&self.handle[0],
            quadtree_handle=&c_quadtree_handle[0]
        )