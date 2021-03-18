module m_cells

    use iso_c_binding

    contains

    subroutine set_2d_computational_nodes_lines(origin, lgrmin, kmax, mmax, nmax, dx,&
        lg, size_i, size_j, nodk, nodm, nodn, quad_idx, bounds, coords, pixel_coords, size_n,&
        area_mask, size_a, size_b, line, n_line_u, n_line_v) bind(c, name="f_set_2d_computational_nodes_lines")
        !!! Entry point for setting nodes and lines and there necessary attributes.
        use m_grid_utils, only : get_lg_corners, get_cell_bbox, get_pix_corners

        integer(kind=c_int), intent(in) :: size_a
        integer(kind=c_int), intent(in) :: size_b
        integer(kind=c_int), intent(in) :: size_i
        integer(kind=c_int), intent(in) :: size_j
        integer(kind=c_int), intent(in) :: size_n
        real(kind=c_double), intent(in) :: origin(2) ! Origin of Quadtree grid
        integer(kind=c_int), intent(in) :: lgrmin ! Number of pixels in cell of smallest refinement level
        integer(kind=c_int), intent(in) :: kmax ! Maximum refinement levels
        integer(kind=c_int), intent(in) :: mmax(kmax) ! X Dimension of each refinement level
        integer(kind=c_int), intent(in) :: nmax(kmax) ! Y Dimension of each refinement level
        real(kind=c_double), intent(in) :: dx(kmax) ! Cell size of each refinement level
        integer(kind=c_int), intent(inout) :: nodk(size_n) ! Array with refinement level of comp node
        integer(kind=c_int), intent(inout) :: nodm(size_n) ! Array with x or m coordinate on its refinement level grid
        integer(kind=c_int), intent(inout) :: nodn(size_n) ! Array with y or n coordinate on its refinement level grid
        integer(kind=c_int), intent(inout) :: lg(size_i,size_j) ! Array with all refinement levels.
        integer(kind=c_int), intent(inout) :: quad_idx(size_i,size_j) ! Array with idx of cell at lg refinement locations
        real(kind=c_double), intent(inout) :: bounds(size_n, 4) ! Bbox of comp cell
        real(kind=c_double), intent(inout) :: coords(size_n, 2) ! Cell center coordinates
        integer(kind=c_int), intent(inout) :: pixel_coords(size_n, 4) ! pixel bbox of comp cell
        integer(kind=c_int), intent(inout) :: area_mask(size_a,size_b) ! Array with active pixels of model.
        integer(kind=c_int), intent(in) :: n_line_u ! Number of active u-dir lines.
        integer(kind=c_int), intent(in) :: n_line_v  ! Number of active v-dir lines.
        integer(kind=c_int), intent(inout) :: line(n_line_u+n_line_v, 2) ! Array with connecting nodes of line.
        integer :: nod
        integer :: k
        integer :: max_pix
        integer :: i0, i1, j0, j1
        integer :: l_u, l_v
        integer :: m, n
        integer :: mn(4)

        max_pix = lgrmin * nmax(1) ! Add 1 because we start indexing at 1 
        write(*,*) '** INFO: Start setting 2D calculation cells.'
        nod = 1
        l_u = 0
        l_v = n_line_u
        line = 0
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
                        call get_pix_corners(k, m, n, lgrmin, i0, i1, j0, j1)
                        pixel_coords(nod, :) = (/ i0 - 1, max_pix - j1, i1 - 1, max_pix - j0 /) ! We inverse the y-axis for pixel_coords to comply with geotiffs in future use.
                        call set_2d_computational_lines(l_u, l_v, k, m, n, mn, lg, lgrmin, area_mask, quad_idx, nod, line)
                        nod = nod + 1
                    else
                        continue
                    endif
                enddo
            enddo
        enddo
        write(*,*) '** INFO: Number of 2D nodes is: ', nod - 1
        write(*,*) '** INFO: Number of 2D lines is: ', l_u + l_v
        write(*,*) '** INFO: Done setting 2D calculation cells.'

    end subroutine set_2d_computational_nodes_lines

    subroutine set_2d_computational_lines(l_u, l_v, k, m, n, mn, lg, lgrmin, area_mask, quad_idx, nod, line)

        use m_grid_utils, only : get_lg_corners, get_pix_corners, crop_pix_coords_to_raster

        integer, intent(inout) :: l_u
        integer, intent(inout) :: l_v
        integer, intent(in) :: k
        integer, intent(in) :: m
        integer, intent(in) :: n
        integer, intent(in) :: mn(4)
        integer, intent(in) :: lg(:,:)
        integer, intent(in) :: lgrmin
        integer, intent(in) :: area_mask(:,:)
        integer, intent(in) :: quad_idx(:,:)
        integer, intent(in), optional :: nod
        integer, intent(inout), optional :: line(:,:)
        integer :: i0,i1,j0,j1,i2,i3,j2,j3
        integer :: neighbour

        !! U direction

        call get_pix_corners(k, m, n, lgrmin, i0, i1, j0, j1, i2, i3, j2, j3)
        call crop_pix_coords_to_raster(area_mask, i0, i1, j0, j1, i2, i3, j2, j3)

        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! U - direction !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!
        if (mn(3) < size(quad_idx, 1)) then
            if(lg(mn(3)+1, mn(2)) == k .and. any(minval(area_mask(i1:i1+1,j0:j1), 1) > 0)) then
                l_u = l_u + 1
                if (present(line)) then
                    neighbour = quad_idx(mn(3)+1, mn(2)) 
                    ! We substract 1 from index to comply with C/python indexing.
                    line(l_u,:) = (/ nod - 1, neighbour - 1/)
                endif
            else
                if (lg(mn(3)+1, mn(2)) == k-1 .and. any(minval(area_mask(i1:i1+1,j0:j2), 1) > 0)) then
                    l_u = l_u + 1
                    if (present(line)) then
                        neighbour = quad_idx(mn(3)+1, mn(2)) 
                        ! We substract 1 from index to comply with C/python indexing.
                        line(l_u,:) = (/ nod - 1, neighbour - 1 /)
                    endif
                endif
                if (lg(mn(3)+1, mn(4)) == k-1 .and. any(minval(area_mask(i1:i1+1,j3:j1), 1) > 0)) then
                    l_u = l_u + 1
                    if (present(line)) then
                        neighbour = quad_idx(mn(3)+1,mn(4))
                        ! We substract 1 from index to comply with C/python indexing.
                        line(l_u,:) = (/ nod - 1, neighbour - 1 /)
                    endif
                endif
            endif
        endif

        if (mn(1) > 1) then
            if(lg(mn(1)-1, mn(2)) == k-1 .and. any(minval(area_mask(i0-1:i0,j0:j2), 1) > 0)) then
                l_u = l_u + 1
                if (present(line)) then
                    neighbour = quad_idx(mn(1)-1, mn(2))
                    ! We substract 1 from index to comply with C/python indexing.
                    line(l_u,:) = (/ neighbour - 1, nod - 1 /)
                endif
            endif
            if(lg(mn(1)-1, mn(4)) == k-1 .and. any(minval(area_mask(i0-1:i0,j3:j1), 1) > 0)) then
                l_u = l_u + 1
                if (present(line)) then
                    neighbour = quad_idx(mn(1)-1, mn(4))
                    ! We substract 1 from index to comply with C/python indexing.
                    line(l_u,:) = (/ neighbour - 1, nod - 1 /)
                endif
            endif
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! U - direction !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! V - direction !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!
        if (mn(4) < size(quad_idx, 2)) then
            if (lg(mn(1), mn(4)+1) == k .and. any(minval(area_mask(i0:i1,j1:j1+1), 2) > 0)) then
                l_v = l_v + 1
                if (present(line)) then
                    neighbour = quad_idx(mn(1), mn(4)+1)
                    ! We substract 1 from index to comply with C/python indexing.
                    line(l_v,:) = (/ nod - 1, neighbour - 1 /)
                endif
            else
                if (lg(mn(1), mn(4)+1) == k-1 .and. any(minval(area_mask(i0:i2,j1:j1+1), 2) > 0)) then
                    l_v = l_v + 1
                    if (present(line)) then
                        neighbour = quad_idx(mn(1), mn(4)+1)
                        ! We substract 1 from index to comply with C/python indexing.
                        line(l_v,:) = (/ nod - 1, neighbour - 1 /)
                    endif
                endif
                if (lg(mn(3), mn(4)+1) == k-1 .and. any(minval(area_mask(i3:i1,j1:j1+1), 2) > 0)) then
                    l_v = l_v + 1
                    if (present(line)) then
                        neighbour = quad_idx(mn(3), mn(4)+1)
                        ! We substract 1 from index to comply with C/python indexing.
                        line(l_v,:) = (/ nod - 1, neighbour - 1 /)
                    endif
                endif
            endif
        endif
            
        if (mn(2) > 1) then
            if(lg(mn(1), mn(2)-1) == k-1 .and. any(minval(area_mask(i0:i2,max(1,j0-1):j0), 2) > 0)) then
                l_v = l_v + 1
                if (present(line)) then
                    neighbour = quad_idx(mn(1), mn(2)-1)
                    ! We substract 1 from index to comply with C/python indexing.
                    line(l_v,:) = (/ neighbour - 1, nod - 1 /)
                endif
            endif
            if(lg(mn(1), mn(2)-1) == k-1 .and. any(minval(area_mask(i3:i1,max(1,j0-1):j0), 2) > 0)) then
                l_v = l_v + 1
                if (present(line)) then
                    neighbour = quad_idx(mn(3), mn(2)-1)
                    ! We substract 1 from index to comply with C/python indexing.
                    line(l_v,:) = (/ neighbour - 1, nod - 1 /)
                endif
            endif
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! V - direction !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine set_2d_computational_lines


end module m_cells