module m_nodes1d

    use parameters, only : NODE_2D,&
                           NODE_2D_GW,&
                           NODE_1D,&
                           NODE_1D_STOR,&
                           NODE_2D_BOUND,&
                           NODE_2D_GW_BOUND,&
                           NODE_1D_BOUND
    use m_nodes
    

    type, extends(ClassNodes) :: ClassNodes1D
        integer, allocatable :: conn_node_id(:)
        integer, allocatable :: conn_node_type(:)
    end type ClassNodes1D

end module m_nodes1d