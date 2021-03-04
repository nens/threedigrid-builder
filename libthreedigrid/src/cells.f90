module m_cells

    use iso_c_binding

    contains

    subroutine set_2d_computational_nodes(origin, kmax, mmax, nmax, dx, size_i, size_j,&
        nodk, nodm, nodn, quad_idx, bounds, coords, size_n) bind(c, name="f_set_2d_computational_nodes")

        use m_grid_utils, only : get_lg_corners, get_cell_bbox

        integer(kind=c_int), intent(in) :: size_i
        integer(kind=c_int), intent(in) :: size_j
        integer(kind=c_int), intent(in) :: size_n
        real(kind=c_double), intent(in) :: origin(2)
        integer(kind=c_int), intent(in) :: kmax
        integer(kind=c_int), intent(in) :: mmax(kmax)
        integer(kind=c_int), intent(in) :: nmax(kmax)
        real(kind=c_double), intent(in) :: dx(kmax)
        integer(kind=c_int), intent(inout) :: nodk(size_n)
        integer(kind=c_int), intent(inout) :: nodm(size_n)
        integer(kind=c_int), intent(inout) :: nodn(size_n)
        integer(kind=c_int), intent(inout) :: quad_idx(size_i,size_j)
        real(kind=c_double), intent(inout) :: bounds(size_n, 4)
        real(kind=c_double), intent(inout) :: coords(size_n, 2)
        integer :: nod
        integer :: k
        integer :: m, n
        integer :: mn(4)

        write(*,*) '** INFO: Start setting 2D calculation cells.'
        nod = 1
        do k=kmax,1,-1
            do m=1,mmax(k)
                do n=1,nmax(k)
                    mn = get_lg_corners(k, m, n)
                    if(all(quad_idx(mn(1):mn(3),mn(2):mn(4)) == nod)) then
                        nodk(nod) = k
                        nodm(nod) = m
                        nodn(nod) = n
                        bounds(nod,:) = get_cell_bbox(origin(1), origin(2), m, n, dx(k))
                        coords(nod, :) = (/ 0.5d0 * (bounds(nod,1) + bounds(nod,3)), 0.5d0 * (bounds(nod,2) + bounds(nod,4)) /)
                        nod = nod + 1
                    else
                        continue
                    endif
                enddo
            enddo
        enddo
        write(*,*) '** INFO: Number of 2D nodes is: ', nod - 1
        write(*,*) '** INFO: Done setting 2D calculation cells.'

    end subroutine set_2d_computational_nodes

    subroutine set_2d_computational_lines(nodk, nodm, nodn, size_m, lgrmin, model_area, size_i, size_j,&
        quad_idx, size_k, size_l, line, size_n) bind(c, name="f_set_2d_computational_lines")

        use m_grid_utils, only : get_lg_corners, get_pix_corners, crop_pix_coords_to_raster

        integer(kind=c_int), intent(in) :: size_m
        integer(kind=c_int), intent(in) :: size_n
        integer(kind=c_int), intent(in) :: size_i
        integer(kind=c_int), intent(in) :: size_j
        integer(kind=c_int), intent(in) :: size_k
        integer(kind=c_int), intent(in) :: size_l
        integer(kind=c_int), intent(in) :: nodk(size_m)
        integer(kind=c_int), intent(in) :: nodm(size_m)
        integer(kind=c_int), intent(in) :: nodn(size_m)
        integer(kind=c_int), intent(in) :: lgrmin
        integer(kind=c_int), intent(in) :: model_area(size_i, size_j)
        integer(kind=c_int), intent(in) :: quad_idx(size_k, size_l)
        integer(kind=c_int), intent(inout) :: line(size_n, 2)
        integer :: nod
        integer :: k, m, n
        integer :: mn(4)
        integer :: i0,i1,j0,j1,i2,i3,j2,j3
        integer :: l
        integer :: neighbour

        l = 0
        line = 0
        !! U direction
        do nod=1,size_m
            k = nodk(nod)
            m = nodm(nod)
            n = nodn(nod)
            mn = get_lg_corners(k, m, n)
            call get_pix_corners(k, m, n, lgrmin, i0, i1, j0, j1, i2, i3, j2, j3)
            call crop_pix_coords_to_raster(model_area, i0, i1, j0, j1, i2, i3, j2, j3)

            !!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!! U - direction !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!
            if (mn(3) < size(quad_idx, 1)) then
                neighbour = quad_idx(mn(3)+1, mn(2)) 
                if (nodk(neighbour) == k.and.any(minval(model_area(i1:i1+1,j0:j1), 1) > 0)) then
                    l = l + 1
                    line(l,:) = (/ nod - 1, neighbour - 1/)
                else
                    !neighbour = quad_idx(m1+1,n0)
                    if (nodk(neighbour) == k-1 .and. any(minval(model_area(i1:i1+1,j0:j2), 1) > 0)) then
                        l = l + 1
                        line(l,:) = (/ nod - 1, neighbour - 1 /)
                    endif
                    neighbour = quad_idx(mn(3)+1,mn(4))
                    if (nodk(neighbour) == k-1 .and. any(minval(model_area(i1:i1+1,j3:j1), 1) > 0)) then
                        l = l + 1
                        line(l,:) = (/ nod - 1, neighbour - 1 /)
                    endif
                endif
            endif

            if (mn(1) > 1) then
                neighbour = quad_idx(mn(1)-1, mn(2))
                if(nodk(neighbour) == k-1 .and. any(minval(model_area(i0-1:i0,j0:j2), 1) > 0)) then
                    l = l + 1
                    line(l,:) = (/ neighbour - 1, nod - 1 /)
                endif
                neighbour = quad_idx(mn(1)-1, mn(4))
                if(nodk(neighbour) == k-1 .and. any(minval(model_area(i0-1:i0,j3:j1), 1) > 0)) then
                    l = l + 1
                    line(l,:) = (/ neighbour - 1, nod - 1 /)
                endif
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!! U - direction !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!


            !!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!! V - direction !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!
            if (mn(4) < size(quad_idx, 2)) then
                neighbour = quad_idx(mn(1), mn(4)+1)
                if (nodk(neighbour) == k .and. any(minval(model_area(i0:i1,j1:j1+1), 2) > 0)) then
                    l = l + 1
                    line(l,:) = (/ nod - 1, neighbour - 1 /)
                else
                    !neighbour = quad_idx(m0, n1+1)
                    if (nodk(neighbour) == k-1 .and. any(minval(model_area(i0:i2,j1:j1+1), 2) > 0)) then
                        l = l + 1
                        line(l,:) = (/ nod - 1, neighbour - 1 /)
                    endif
                    neighbour = quad_idx(mn(3), mn(4)+1)
                    if (nodk(neighbour) == k-1 .and. any(minval(model_area(i3:i1,j1:j1+1), 2) > 0)) then
                        l = l + 1
                        line(l,:) = (/ nod - 1, neighbour - 1 /)
                    endif
                endif
            endif
            
            if (mn(2) > 1) then
                neighbour = quad_idx(mn(1), mn(2)-1)
                if(nodk(neighbour) == k-1 .and. any(minval(model_area(i0:i2,max(1,j0-1):j0), 2) > 0)) then
                    l = l + 1
                    line(l,:) = (/ neighbour - 1, nod - 1 /)
                endif
                neighbour = quad_idx(mn(3), mn(2)-1)
                if(nodk(neighbour) == k-1 .and. any(minval(model_area(i3:i1,max(1,j0-1):j0), 2) > 0)) then
                    l = l + 1
                    line(l,:) = (/ neighbour - 1, nod - 1 /)
                endif
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!! V - direction !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!
        enddo
        write(*,*) '** INFO: Number of 2D Surface flow lines is: ', l

    end subroutine set_2d_computational_lines


end module m_cells