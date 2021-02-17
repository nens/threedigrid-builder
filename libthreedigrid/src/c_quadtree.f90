module m_c_quadtree

    use iso_c_binding

    implicit none

    contains

    subroutine c_set_refinement(refine_id, refine_geom, n0, n1, refine_level,&
        refine_type, origin, bbox, mmax, nmax, dx, j0, lg, i0, i1)&
        bind(c, name="f_set_refinement")
        
        use m_quadtree, only : set_refinement

        integer(kind=c_int), intent(in) :: n0
        integer(kind=c_int), intent(in) :: n1
        integer(kind=c_int), intent(in) :: j0
        integer(kind=c_int), intent(in) :: i0
        integer(kind=c_int), intent(in) :: i1
        real(kind=c_double), intent(in) :: refine_geom(n0,n1)
        integer(kind=c_int), intent(in) :: refine_id
        integer(kind=c_int), intent(in) :: refine_level
        integer(kind=c_int), intent(in) :: refine_type
        real(kind=c_double), intent(in) :: origin(2)
        real(kind=c_double), intent(in) :: bbox(4)
        integer(kind=c_int), intent(in) :: mmax(j0)
        integer(kind=c_int), intent(in) :: nmax(j0)
        real(kind=c_double), intent(in) :: dx(j0)
        integer(kind=c_int), intent(inout) :: lg(i0,i1)

        call set_refinement(refine_id=refine_id,&
                            refine_geom=refine_geom,&
                            refine_level=refine_level,&
                            refine_type=refine_type,&
                            origin=origin,&
                            bbox=bbox,&
                            mmax=mmax,&
                            nmax=nmax,&
                            dx=dx,&
                            lg=lg)

    end subroutine c_set_refinement

    subroutine c_make_quadtree(kmax, mmax, nmax, lg, n0, i0, i1) bind(c, name="make_quadtree")

        use m_quadtree, only : make_quadtree

        integer(kind=c_int), intent(in) :: n0
        integer(kind=c_int), intent(in) :: i0
        integer(kind=c_int), intent(in) :: i1
        integer(kind=c_int), intent(in) :: kmax
        integer(kind=c_int), intent(in) :: mmax(n0)
        integer(kind=c_int), intent(in) :: nmax(n0)
        integer(kind=c_int), intent(inout) :: lg(i0,i1)

        call make_quadtree(kmax=kmax, mmax=mmax, nmax=nmax, lg=lg)

    end subroutine c_make_quadtree

    ! subroutine c_set_active_2d_comp_cells(handle, model_area, n0, n1) bind(c, name="set_active_2d_comp_cells")
! 
        ! use m_quadtree, only : QuadTreeFortran, set_active_2d_comp_cells
! 
        ! type quadtreefortran_ptr_type
            ! type(QuadTreeFortran), pointer :: p => NULL()
        ! end type quadtreefortran_ptr_type
        ! type(quadtreefortran_ptr_type) :: handle_ptr
        ! integer(kind=c_int), intent(inout) :: handle(2)
        ! integer(kind=c_int), intent(in) :: n0
        ! integer(kind=c_int), intent(in) :: n1
        ! integer(kind=c_int), intent(in) :: model_area(n0,n1)
! 
        ! handle_ptr = transfer(handle, handle_ptr)
        ! call set_active_2d_comp_cells(self=handle_ptr%p, model_area=model_area)
! 
    ! end subroutine c_set_active_2d_comp_cells

end module m_c_quadtree