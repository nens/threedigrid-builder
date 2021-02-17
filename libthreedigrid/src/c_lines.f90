! module m_c_lines

!     use iso_c_binding

!     implicit none

!     contains

!     subroutine c_init_lines(handle, nodes_handle, model_area, n0, n1) bind(c, name="init_lines")

!         use m_lines, only : ClassLines, init_lines
!         use m_nodes, only : ClassNodes

!         type classlines_ptr_type
!             type(ClassLines), pointer :: p => NULL()
!         end type classlines_ptr_type
!         type(classlines_ptr_type) :: handle_ptr
!         integer(kind=c_int), intent(out) :: handle(2)
!         type classnodes_ptr_type
!             type(ClassNodes), pointer :: p => NULL()
!         end type classnodes_ptr_type
!         type(classnodes_ptr_type) :: nodes_handle_ptr
!         integer(kind=c_int), intent(in) :: nodes_handle(2)
!         integer(kind=c_int), intent(in) :: n0
!         integer(kind=c_int), intent(in) :: n1
!         integer(kind=c_int), intent(in) :: model_area(n0,n1)

!         nodes_handle_ptr = transfer(nodes_handle, nodes_handle_ptr)
!         allocate(handle_ptr%p)
!         call init_lines(self=handle_ptr%p, nodes=nodes_handle_ptr%p, model_area=model_area)
!         handle = transfer(handle_ptr, handle)
        
!     end subroutine c_init_lines


!     function c_get_lintot(handle) result(lintot) bind(c, name="get_lintot")

!         use m_lines, only : ClassLines, get_lintot

!         type classlines_ptr_type
!             type(ClassLines), pointer :: p => NULL()
!         end type classlines_ptr_type
!         type(classlines_ptr_type) :: handle_ptr
!         integer(kind=c_int), intent(inout) :: handle(2)
!         integer :: lintot

!         handle_ptr = transfer(handle, handle_ptr)
!         lintot = get_lintot(self=handle_ptr%p)

!     end function c_get_lintot


!     function c_get_id(handle, l0, l1) result(id_cnt) bind(c, name="get_line_id")

!         use m_lines, only : ClassLines, get_id

!         type classlines_ptr_type
!             type(ClassLines), pointer :: p => NULL()
!         end type classlines_ptr_type
!         type(classlines_ptr_type) :: handle_ptr
!         integer(kind=c_int), intent(inout) :: handle(2)
!         integer(kind=c_int), intent(in) :: l0
!         integer(kind=c_int), intent(in) :: l1
!         integer, pointer :: ids(:)
!         type(C_PTR) :: id_cnt

!         handle_ptr = transfer(handle, handle_ptr)
!         ids => get_id(self=handle_ptr%p, l0=l0, l1=l1)
!         id_cnt = c_loc(ids)
!         nullify(ids)

!     end function c_get_id


!     function c_get_line(handle, l0, l1) result(line_cnt) bind(c, name="get_line")

!         use m_lines, only : ClassLines, get_line

!         type classlines_ptr_type
!             type(ClassLines), pointer :: p => NULL()
!         end type classlines_ptr_type
!         type(classlines_ptr_type) :: handle_ptr
!         integer(kind=c_int), intent(inout) :: handle(2)
!         integer(kind=c_int), intent(in) :: l0
!         integer(kind=c_int), intent(in) :: l1
!         integer, pointer :: line(:,:)
!         type(C_PTR) :: line_cnt

!         handle_ptr = transfer(handle, handle_ptr)
!         line => get_line(self=handle_ptr%p, l0=l0, l1=l1)
!         line_cnt = c_loc(line)
!         nullify(line)

!     end function c_get_line

! end module m_c_lines