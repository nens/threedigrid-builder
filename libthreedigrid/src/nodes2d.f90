module m_nodes2d

    use m_nodes
    use parameters, only : NODE_2D,&
                           NODE_2D_GW,&
                           NODE_1D,&
                           NODE_1D_STOR,&
                           NODE_2D_BOUND,&
                           NODE_2D_GW_BOUND,&
                           NODE_1D_BOUND

    type, extends(ClassNodes) :: ClassNodes2D
    end type ClassNodes2D

    contains



end module m_nodes2d