module m_cells

    use iso_c_binding

    contains

    subroutine set_2d_computational_nodes(origin, kmax, mmax, nmax, dx, lg, size_i, size_j,&
        nodk, nodm, nodn, quad_nod, bounds, coords, size_n) bind(c, name="f_set_2d_computational_nodes")

        use m_grid_utils, only : get_lg_corners, get_cell_bbox

        integer(kind=c_int), intent(in) :: size_i
        integer(kind=c_int), intent(in) :: size_j
        integer(kind=c_int), intent(in) :: size_n
        real(kind=c_double), intent(in) :: origin(2)
        integer(kind=c_int), intent(in) :: kmax
        integer(kind=c_int), intent(in) :: mmax(kmax)
        integer(kind=c_int), intent(in) :: nmax(kmax)
        real(kind=c_double), intent(in) :: dx(kmax)
        integer(kind=c_int), intent(in) :: lg(size_i,size_j)
        integer(kind=c_int), intent(inout) :: nodk(size_n)
        integer(kind=c_int), intent(inout) :: nodm(size_n)
        integer(kind=c_int), intent(inout) :: nodn(size_n)
        integer(kind=c_int), intent(inout) :: quad_nod(size_i,size_j)
        real(kind=c_double), intent(inout) :: bounds(size_n, 4)
        real(kind=c_double), intent(inout) :: coords(size_n, 2)
        integer :: nod
        integer :: k
        integer :: m, n
        integer :: m0, m1, n0, n1

        write(*,*) '** INFO: Start setting 2D calculation cells.'
        nod = 0
        do k=kmax,1,-1
            do m=1,mmax(k)
                do n=1,nmax(k)
                    call get_lg_corners(k, m, n, m0, m1, n0, n1)
                    if(all(lg(m0:m1,n0:n1) == k)) then
                        nod = nod + 1
                        nodk(nod) = k
                        nodm(nod) = m
                        nodn(nod) = n
                        quad_nod(m0:m1,n0:n1) = nod
                        bounds(nod,:) = get_cell_bbox(origin(1), origin(2), m, n, dx(k))
                        coords(nod, :) = (/ 0.5d0 * (bounds(nod,1) + bounds(nod,3)), 0.5d0 * (bounds(nod,2) + bounds(nod,4)) /)
                    else
                        continue
                    endif
                enddo
            enddo
        enddo
        write(*,*) '** INFO: Number of 2D nodes is: ', nod
        write(*,*) '** INFO: Done setting 2D calculation cells.'

    end subroutine set_2d_computational_nodes

end module m_cells