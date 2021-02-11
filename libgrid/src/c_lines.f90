module m_c_lines

    use iso_c_binding

    implicit none

    contains

    subroutine c_init_lines(handle, nodes_handle, model_area, n0, n1) bind(c, name="init_lines")

        use m_lines, only : ClassLines, init_lines
        use m_nodes, only : ClassNodes

        type classlines_ptr_type
            type(ClassLines), pointer :: p => NULL()
        end type classlines_ptr_type
        type(classlines_ptr_type) :: handle_ptr
        integer(kind=c_int), intent(out) :: handle(2)
        type classnodes_ptr_type
            type(ClassNodes), pointer :: p => NULL()
        end type classnodes_ptr_type
        type(classnodes_ptr_type) :: nodes_handle_ptr
        integer(kind=c_int), intent(in) :: nodes_handle(2)
        integer(kind=c_int), intent(in) :: n0
        integer(kind=c_int), intent(in) :: n1
        integer(kind=c_int), intent(in) :: model_area(n0,n1)

        nodes_handle_ptr = transfer(nodes_handle, nodes_handle_ptr)
        allocate(handle_ptr%p)
        call init_lines(self=handle_ptr%p, nodes=nodes_handle_ptr%p, model_area=model_area)
        handle = transfer(handle_ptr, handle)
        
    end subroutine c_init_lines

end module m_c_lines