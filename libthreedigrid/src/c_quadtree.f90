module m_c_quadtree

    use iso_c_binding

    implicit none

    contains

    subroutine c_init_quadtree(handle, x0p, y0p, dxp, imax, jmax, lgrmin, kmax)&
        bind(c, name="init_quadtree")

        use m_quadtree, only : QuadTreeFortran, init_quadtree

        type quadtreefortran_ptr_type
            type(QuadTreeFortran), pointer :: p => NULL()
        end type quadtreefortran_ptr_type
        type(quadtreefortran_ptr_type) :: handle_ptr
        integer(kind=c_int), intent(out) :: handle(2)
        real(kind=c_double), intent(in) :: x0p
        real(kind=c_double), intent(in) :: y0p
        real(kind=c_double), intent(in) :: dxp
        integer(kind=c_int), intent(in) :: imax
        integer(kind=c_int), intent(in) :: jmax
        integer(kind=c_int), intent(in) :: lgrmin
        integer(kind=c_int), intent(in) :: kmax

        write(*,*) 'OH HELLOOOOO !!!!!'
        allocate(handle_ptr%p)
        call init_quadtree(self=handle_ptr%p, x0p=x0p, y0p=y0p, dxp=dxp, imax=imax, jmax=jmax, lgrmin=lgrmin, kmax=kmax)
        handle = transfer(handle_ptr, handle)
        
    end subroutine c_init_quadtree

    function c_set_refinement(handle, refine_id, refine_geom, n0, n1, refine_level, refine_type) result(status)&
        bind(c, name="set_refinement")

        use m_quadtree, only : QuadTreeFortran, set_refinement

        type quadtreefortran_ptr_type
            type(QuadTreeFortran), pointer :: p => NULL()
        end type quadtreefortran_ptr_type
        type(quadtreefortran_ptr_type) :: handle_ptr
        integer(kind=c_int), intent(inout) :: handle(2)
        integer(kind=c_int), intent(in) :: n0
        integer(kind=c_int), intent(in) :: n1
        real(kind=c_double), intent(in) :: refine_geom(n0,n1)
        integer(kind=c_int), intent(in) :: refine_id
        integer(kind=c_int), intent(in) :: refine_level
        integer(kind=c_int), intent(in) :: refine_type
        integer(kind=c_int) :: status

        handle_ptr = transfer(handle, handle_ptr)
        call set_refinement(self=handle_ptr%p,&
                            refine_id=refine_id,&
                            refine_geom=refine_geom,&
                            refine_level=refine_level,&
                            refine_type=refine_type,&
                            status=status)

    end function c_set_refinement

    subroutine c_make_quadtree(handle) bind(c, name="make_quadtree")

        use m_quadtree, only : QuadTreeFortran, make_quadtree

        type quadtreefortran_ptr_type
            type(QuadTreeFortran), pointer :: p => NULL()
        end type quadtreefortran_ptr_type
        type(quadtreefortran_ptr_type) :: handle_ptr
        integer(kind=c_int), intent(out) :: handle(2)

        handle_ptr = transfer(handle, handle_ptr)
        call make_quadtree(self=handle_ptr%p)

    end subroutine c_make_quadtree

    subroutine c_set_active_2d_comp_cells(handle, model_area, n0, n1) bind(c, name="set_active_2d_comp_cells")

        use m_quadtree, only : QuadTreeFortran, set_active_2d_comp_cells

        type quadtreefortran_ptr_type
            type(QuadTreeFortran), pointer :: p => NULL()
        end type quadtreefortran_ptr_type
        type(quadtreefortran_ptr_type) :: handle_ptr
        integer(kind=c_int), intent(inout) :: handle(2)
        integer(kind=c_int), intent(in) :: n0
        integer(kind=c_int), intent(in) :: n1
        integer(kind=c_int), intent(in) :: model_area(n0,n1)

        handle_ptr = transfer(handle, handle_ptr)
        call set_active_2d_comp_cells(self=handle_ptr%p, model_area=model_area)

    end subroutine c_set_active_2d_comp_cells

end module m_c_quadtree