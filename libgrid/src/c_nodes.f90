module m_c_nodes

    use iso_c_binding

    implicit none

    contains

    subroutine c_init_nodes(handle, quadtree_handle) bind(c, name="init_nodes")

        use m_quadtree, only : QuadTreeFortran
        use m_nodes, only : ClassNodes, init_nodes

        type classnodes_ptr_type
            type(ClassNodes), pointer :: p => NULL()
        end type classnodes_ptr_type
        type(classnodes_ptr_type) :: handle_ptr
        integer(kind=c_int), intent(out) :: handle(2)
        type quadtreefortran_ptr_type
            type(QuadTreeFortran), pointer :: p => NULL()
        end type quadtreefortran_ptr_type
        type(quadtreefortran_ptr_type) :: quadtree_handle_ptr
        integer(kind=c_int), intent(in) :: quadtree_handle(2)

        quadtree_handle_ptr = transfer(quadtree_handle, quadtree_handle_ptr)
        allocate(handle_ptr%p)
        call init_nodes(self=handle_ptr%p, quadtree_grid=quadtree_handle_ptr%p)
        handle = transfer(handle_ptr, handle)
        
    end subroutine c_init_nodes


    function c_get_node_bnds(handle, n0, n1) result(node_bnds_cnt) bind(c, name="get_node_bnds")

        use m_nodes, only : ClassNodes, get_node_bnds

        type classnodes_ptr_type
            type(ClassNodes), pointer :: p => NULL()
        end type classnodes_ptr_type
        type(classnodes_ptr_type) :: handle_ptr
        integer(kind=c_int), intent(inout) :: handle(2)
        integer(kind=c_int), intent(in) :: n0
        integer(kind=c_int), intent(in) :: n1
        double precision, target :: node_bnds(n0:n1, 4)
        type(C_PTR) :: node_bnds_cnt

        handle_ptr = transfer(handle, handle_ptr)
        node_bnds = get_node_bnds(self=handle_ptr%p, n_start=n0, n_end=n1)
        node_bnds_cnt = c_loc(node_bnds)
        write(*,*) node_bnds(:,:)

    end function c_get_node_bnds


end module m_c_nodes