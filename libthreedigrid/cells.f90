module m_cells

    implicit none

    contains

    subroutine set_2d_computational_nodes_lines(origin, lgrmin, kmax, mmax, nmax, dx,&
        lg, nodk, nodm, nodn, quad_idx, bounds, coords, pixel_coords,&
        area_mask, line, cross_pix_coords, n_line_u, n_line_v)
        !!! Entry point for setting nodes and lines and there necessary attributes.
        use m_grid_utils, only : get_lg_corners, get_cell_bbox, get_pix_corners, pad_area_mask

        double precision, intent(in) :: origin(2) ! Origin of Quadtree grid
        integer, intent(in) :: lgrmin ! Number of pixels in cell of smallest refinement level
        integer, intent(in) :: kmax ! Maximum refinement levels
        integer, intent(in) :: mmax(:) ! X Dimension of each refinement level
        integer, intent(in) :: nmax(:) ! Y Dimension of each refinement level
        double precision, intent(in) :: dx(:) ! Cell size of each refinement level
        integer, intent(inout) :: nodk(:) ! Array with refinement level of comp node
        integer, intent(inout) :: nodm(:) ! Array with x or m coordinate on its refinement level grid
        integer, intent(inout) :: nodn(:) ! Array with y or n coordinate on its refinement level grid
        integer, intent(inout) :: lg(:, :) ! Array with all refinement levels.
        integer, intent(inout) :: quad_idx(:, :) ! Array with idx of cell at lg refinement locations
        double precision, intent(inout) :: bounds(:, :) ! Bbox of comp cell
        double precision, intent(inout) :: coords(:, :) ! Cell center coordinates
        integer, intent(inout) :: pixel_coords(:, :) ! pixel bbox of comp cell
        integer*1, intent(inout) :: area_mask(:, :) ! Array with active pixels of model.
        integer, intent(inout) :: line(:, :) ! Array with connecting nodes of line.
        integer, intent(inout) :: cross_pix_coords(:, :) ! Array pixel indices of line interface
        integer, intent(in) :: n_line_u ! Number of active u-dir lines.
        integer, intent(in) :: n_line_v  ! Number of active v-dir lines.
        integer*1, allocatable :: area_mask_padded(:, :)
        integer :: nod
        integer :: k
        integer :: i0, i1, j0, j1
        integer :: l_u, l_v
        integer :: m, n
        integer :: mn(4)
        logical :: use_2d_flow

        write(*,*) '** INFO: Start setting 2D calculation cells.'
        nod = 1
        use_2d_flow = ((n_line_u > 0) .or. (n_line_v > 0))
        l_u = 0
        l_v = n_line_u
        line = 0
        call get_pix_corners(kmax, mmax(kmax), nmax(kmax), lgrmin, i0, i1, j0, j1)
        if (use_2d_flow) then
            area_mask_padded = pad_area_mask(area_mask, i0, i1, j0, j1) 
        endif
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
                        ! We inverse the y-axis for pixel_coords to comply with geotiffs in future use.
                        ! And do some index fiddling because python starts indexing at 0 and has open end indexing.
                        pixel_coords(nod, :) = (/ i0 - 1, j0 - 1, i1, j1 /)
                        if (use_2d_flow) then
                            call set_2d_computational_lines(&
                                l_u, l_v, k, m, n, mn, lg, lgrmin, area_mask_padded, quad_idx, nod, line, cross_pix_coords&
                            )
                        endif
                        nod = nod + 1
                    else
                        continue
                    endif
                enddo
            enddo
        enddo
        if (use_2d_flow) then
            deallocate(area_mask_padded)
        endif
        write(*,*) '** INFO: Number of 2D nodes is: ', nod - 1
        write(*,*) '** INFO: Number of 2D lines is: ', l_u + l_v
        write(*,*) '** INFO: Done setting 2D calculation cells.'

    end subroutine set_2d_computational_nodes_lines

    subroutine set_2d_computational_lines(l_u, l_v, k, m, n, mn, lg, lgrmin, area_mask, quad_idx, nod, line, cross_pix_coords)

        use m_grid_utils, only : get_lg_corners, get_pix_corners, crop_pix_coords_to_raster

        integer, intent(inout) :: l_u
        integer, intent(inout) :: l_v
        integer, intent(in) :: k
        integer, intent(in) :: m
        integer, intent(in) :: n
        integer, intent(in) :: mn(4)
        integer, intent(in) :: lg(:,:)
        integer, intent(in) :: lgrmin
        integer*1, intent(in) :: area_mask(:,:)
        integer, intent(in) :: quad_idx(:,:)
        integer, intent(in), optional :: nod
        integer, intent(inout), optional :: line(:,:)
        integer, intent(inout), optional :: cross_pix_coords(:,:)
        integer :: i0,i1,j0,j1,i2,i3,j2,j3
        integer :: neighbour

        !! U direction

        call get_pix_corners(k, m, n, lgrmin, i0, i1, j0, j1, i2, i3, j2, j3)

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
                    cross_pix_coords(l_u,:) = (/ i1, j0 - 1, i1, j1 /) ! Python indexing so -1 at start slice
                endif
            else
                if (lg(mn(3)+1, mn(2)) == k-1 .and. any(minval(area_mask(i1:i1+1,j0:j2), 1) > 0)) then
                    l_u = l_u + 1
                    if (present(line)) then
                        neighbour = quad_idx(mn(3)+1, mn(2)) 
                        ! We substract 1 from index to comply with C/python indexing.
                        line(l_u,:) = (/ nod - 1, neighbour - 1 /)
                        cross_pix_coords(l_u,:) = (/ i1, j0 - 1, i1, j2 /) ! Python indexing so -1 at start slice
                    endif
                endif
                if (lg(mn(3)+1, mn(4)) == k-1 .and. any(minval(area_mask(i1:i1+1,j3:j1), 1) > 0)) then
                    l_u = l_u + 1
                    if (present(line)) then
                        neighbour = quad_idx(mn(3)+1,mn(4))
                        ! We substract 1 from index to comply with C/python indexing.
                        line(l_u,:) = (/ nod - 1, neighbour - 1 /)
                        cross_pix_coords(l_u,:) = (/ i1, j3 - 1, i1, j1 /) ! Python indexing so -1 at start slice
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
                    cross_pix_coords(l_u,:) = (/ i0 - 1, j0 - 1, i0 - 1, j2 /)
                endif
            endif
            if(lg(mn(1)-1, mn(4)) == k-1 .and. any(minval(area_mask(i0-1:i0,j3:j1), 1) > 0)) then
                l_u = l_u + 1
                if (present(line)) then
                    neighbour = quad_idx(mn(1)-1, mn(4))
                    ! We substract 1 from index to comply with C/python indexing.
                    line(l_u,:) = (/ neighbour - 1, nod - 1 /)
                    cross_pix_coords(l_u,:) = (/ i0 - 1, j3 - 1, i0 - 1, j1 /)
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
                    cross_pix_coords(l_v,:) = (/ i0 - 1, j1, i1, j1 /) ! Python indexing so -1 at start slice
                endif
            else
                if (lg(mn(1), mn(4)+1) == k-1 .and. any(minval(area_mask(i0:i2,j1:j1+1), 2) > 0)) then
                    l_v = l_v + 1
                    if (present(line)) then
                        neighbour = quad_idx(mn(1), mn(4)+1)
                        ! We substract 1 from index to comply with C/python indexing.
                        line(l_v,:) = (/ nod - 1, neighbour - 1 /)
                        cross_pix_coords(l_v,:) = (/ i0 - 1, j1, i2, j1 /) ! Python indexing so -1 at start slice
                    endif
                endif
                if (lg(mn(3), mn(4)+1) == k-1 .and. any(minval(area_mask(i3:i1,j1:j1+1), 2) > 0)) then
                    l_v = l_v + 1
                    if (present(line)) then
                        neighbour = quad_idx(mn(3), mn(4)+1)
                        ! We substract 1 from index to comply with C/python indexing.
                        line(l_v,:) = (/ nod - 1, neighbour - 1 /)
                        cross_pix_coords(l_v,:) = (/ i3 - 1, j1, i1, j1 /) ! Python indexing so -1 at start slice
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
                    cross_pix_coords(l_v,:) = (/ i0 - 1, j0 - 1, i2, j0 - 1/) ! Python indexing so -1 at start slice
                endif
            endif
            if(lg(mn(3), mn(2)-1) == k-1 .and. any(minval(area_mask(i3:i1,max(1,j0-1):j0), 2) > 0)) then
                l_v = l_v + 1
                if (present(line)) then
                    neighbour = quad_idx(mn(3), mn(2)-1)
                    ! We substract 1 from index to comply with C/python indexing.
                    line(l_v,:) = (/ neighbour - 1, nod - 1 /)
                    cross_pix_coords(l_v,:) = (/ i3 - 1, j0 - 1, i1, j0 - 1 /) ! Python indexing so -1 at start slice
                endif
            endif
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! V - direction !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine set_2d_computational_lines


end module m_cells