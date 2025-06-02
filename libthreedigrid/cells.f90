module m_cells

    use m_log

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


        call print_unix('** INFO: Start setting 2D calculation cells.')

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

        call print_unix('** INFO: Number of 2D nodes is: ', nod - 1)
        call print_unix('** INFO: Number of 2D lines is: ', l_u + (l_v - l_u))
        call print_unix('** INFO: Done setting 2D calculation cells.')

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


    subroutine set_quarter_admin(nodk, nodm, nodn, line, kcu, quarter_line, quarter_neighbour, liutot, livtot, n2dobc)

        use parameters, only : LINE_2D_BOUNDARY_EAST, LINE_2D_BOUNDARY_WEST, LINE_2D_BOUNDARY_SOUTH, LINE_2D_BOUNDARY_NORTH

        integer, intent(in) :: nodk(:)
        integer, intent(in) :: nodm(:)
        integer, intent(in) :: nodn(:)
        integer, intent(in) :: line(:, :)
        integer, intent(in) :: kcu(:)
        integer, intent(inout) :: quarter_line(:, :)
        integer, intent(inout) :: quarter_neighbour(:, :)
        integer, intent(in) :: liutot
        integer, intent(in) :: livtot
        integer, intent(in) :: n2dobc
        integer :: l
        integer :: quarter
        integer :: nodd
        integer :: nodu
        integer :: kd
        integer :: ku
        integer :: md
        integer :: mu
        integer :: nd
        integer :: nu
        integer :: nb


        !!!!!!!!!!!!!!!!
        !!  3  !!  4  !!
        !!!!!!!!!!!!!!!!
        !!  1  !!  2  !!
        !!!!!!!!!!!!!!!!
        ! Here we try to find all horizontal flow lines associated with a cell quarter. 
        ! For two adjacent cells of samen size both neighbouring quadrants have same flow line. 
        ! When adjacent cells differ in size due to refinement, each quadrant is associated with its own line for the larger cell. The smaller cell will have the same flow line for both quadrants. Hence the "if" structure.
        ! Same logic applies to finding neighbouring nodes for quarters.
        ! The cell quadrant index is illustrated above.

        do l=1, liutot
            nodd = line(l, 1) + 1 ! This is a 0-based python index, therefore +1
            nodu = line(l, 2) + 1
            kd = nodk(nodd)
            ku = nodk(nodu)
            md = nodm(nodd) 
            mu = nodm(nodu)
            nd = nodn(nodd)
            nu = nodn(nodu)
            if (ku==kd) then
                quarter_line(get_quarter_idx(nodd, 2), 1) = l - 1 ! This needs to be an index in Python again.
                quarter_line(get_quarter_idx(nodd, 4), 1) = l - 1
                quarter_line(get_quarter_idx(nodu, 1), 1) = l - 1
                quarter_line(get_quarter_idx(nodu, 3), 1) = l - 1
                quarter_neighbour(get_quarter_idx(nodd, 2), 1) = nodu - 1
                quarter_neighbour(get_quarter_idx(nodd, 4), 1) = nodu - 1
                quarter_neighbour(get_quarter_idx(nodu, 1), 1) = nodd - 1
                quarter_neighbour(get_quarter_idx(nodu, 3), 1) = nodd - 1
            elseif (ku == kd + 1) then
                if (nd == 2 * nu - 1) then
                    quarter_line(get_quarter_idx(nodd, 2), 1) = l - 1
                    quarter_line(get_quarter_idx(nodd, 4), 1) = l - 1
                    quarter_line(get_quarter_idx(nodu, 1), 1) = l - 1
                    quarter_neighbour(get_quarter_idx(nodd, 2), 1) = nodu - 1
                    quarter_neighbour(get_quarter_idx(nodd, 4), 1) = nodu - 1
                    quarter_neighbour(get_quarter_idx(nodu, 1), 1) = nodd - 1
                endif
                if (nd == 2 * nu) then
                    quarter_line(get_quarter_idx(nodd, 2), 1) = l - 1
                    quarter_line(get_quarter_idx(nodd, 4), 1) = l - 1
                    quarter_line(get_quarter_idx(nodu, 3), 1) = l - 1
                    quarter_neighbour(get_quarter_idx(nodd, 2), 1) = nodu - 1
                    quarter_neighbour(get_quarter_idx(nodd, 4), 1) = nodu - 1
                    quarter_neighbour(get_quarter_idx(nodu, 3), 1) = nodd - 1
                endif
            else if (kd == ku + 1) then
                if (nu == 2 * nd - 1) then
                    quarter_line(get_quarter_idx(nodd, 2), 1) = l - 1
                    quarter_line(get_quarter_idx(nodu, 1), 1) = l - 1
                    quarter_line(get_quarter_idx(nodu, 3), 1) = l - 1
                    quarter_neighbour(get_quarter_idx(nodd, 2), 1) = nodu - 1
                    quarter_neighbour(get_quarter_idx(nodu, 1), 1) = nodd - 1
                    quarter_neighbour(get_quarter_idx(nodu, 3), 1) = nodd - 1
                elseif (nu == 2 * nd) then
                    quarter_line(get_quarter_idx(nodd, 4), 1) = l - 1
                    quarter_line(get_quarter_idx(nodu, 1), 1) = l - 1
                    quarter_line(get_quarter_idx(nodu, 3), 1) = l - 1
                    quarter_neighbour(get_quarter_idx(nodd, 4), 1) = nodu - 1
                    quarter_neighbour(get_quarter_idx(nodu, 1), 1) = nodd - 1
                    quarter_neighbour(get_quarter_idx(nodu, 3), 1) = nodd - 1
                endif
            endif
        enddo

        ! Here we try to find all vertical flow lines associated with a cell quarter. 
        ! For two adjacent cells of samen size both neighbouring quadrants have same flow line. 
        ! When adjacent cells differ in size due to refinement, each quadrant is associated with its own line for the larger cell. The smaller cell will have the same flow line for both quadrants. Hence the "if" structure.
        ! Same logic applies to finding neighbouring nodes for quarters.
        ! The cell quadrant index is illustrated above.
        do l=liutot + 1, liutot+livtot
            nodd = line(l, 1) + 1
            nodu = line(l, 2) + 1 
            kd = nodk(nodd)
            ku = nodk(nodu)
            md = nodm(nodd) 
            mu = nodm(nodu)
            nd = nodn(nodd)
            nu = nodn(nodu)
            if (ku==kd) then
                quarter_line(get_quarter_idx(nodd, 3), 2) = l - 1
                quarter_line(get_quarter_idx(nodd, 4), 2) = l - 1
                quarter_line(get_quarter_idx(nodu, 1), 2) = l - 1
                quarter_line(get_quarter_idx(nodu, 2), 2) = l - 1
                quarter_neighbour(get_quarter_idx(nodd, 3), 2) = nodu - 1
                quarter_neighbour(get_quarter_idx(nodd, 4), 2) = nodu - 1
                quarter_neighbour(get_quarter_idx(nodu, 1), 2) = nodd - 1
                quarter_neighbour(get_quarter_idx(nodu, 2), 2) = nodd - 1
            elseif (ku == kd + 1) then
                if (md == 2 * mu - 1) then
                    quarter_line(get_quarter_idx(nodd, 3), 2) = l - 1
                    quarter_line(get_quarter_idx(nodd, 4), 2) = l - 1
                    quarter_line(get_quarter_idx(nodu, 1), 2) = l - 1
                    quarter_neighbour(get_quarter_idx(nodd, 3), 2) = nodu - 1
                    quarter_neighbour(get_quarter_idx(nodd, 4), 2) = nodu - 1
                    quarter_neighbour(get_quarter_idx(nodu, 1), 2) = nodd - 1
                endif
                if (md == 2 * mu) then
                    quarter_line(get_quarter_idx(nodd, 3), 2) = l - 1
                    quarter_line(get_quarter_idx(nodd, 4), 2) = l - 1
                    quarter_line(get_quarter_idx(nodu, 2), 2) = l - 1
                    quarter_neighbour(get_quarter_idx(nodd, 3), 2) = nodu - 1
                    quarter_neighbour(get_quarter_idx(nodd, 4), 2) = nodu - 1
                    quarter_neighbour(get_quarter_idx(nodu, 2), 2) = nodd - 1
                endif
            else if (kd == ku + 1) then
                if (mu == 2 * md - 1) then
                    quarter_line(get_quarter_idx(nodd, 3), 2) = l - 1
                    quarter_line(get_quarter_idx(nodu, 1), 2) = l - 1
                    quarter_line(get_quarter_idx(nodu, 2), 2) = l - 1
                    quarter_neighbour(get_quarter_idx(nodd, 3), 2) = nodu - 1
                    quarter_neighbour(get_quarter_idx(nodu, 1), 2) = nodd - 1
                    quarter_neighbour(get_quarter_idx(nodu, 2), 2) = nodd - 1
                elseif (mu == 2 * md) then
                    quarter_line(get_quarter_idx(nodd, 4), 2) = l - 1
                    quarter_line(get_quarter_idx(nodu, 1), 2) = l - 1
                    quarter_line(get_quarter_idx(nodu, 2), 2) = l - 1
                    quarter_neighbour(get_quarter_idx(nodd, 4), 2) = nodu - 1
                    quarter_neighbour(get_quarter_idx(nodu, 1), 2) = nodd - 1
                    quarter_neighbour(get_quarter_idx(nodu, 2), 2) = nodd - 1
                endif
            endif
        enddo

        do nb=1,n2dobc  !nodobc
            l = liutot + livtot + nb
            nodd = line(l, 1) + 1
            nodu = line(l, 2) + 1
            select case(kcu(l))
            case(LINE_2D_BOUNDARY_WEST)
                quarter_line(get_quarter_idx(nodu, 1), 1) = l - 1
                quarter_line(get_quarter_idx(nodu, 3), 1) = l
                quarter_neighbour(get_quarter_idx(nodu, 1), 1) = nodd - 1
                quarter_neighbour(get_quarter_idx(nodu, 3), 1) = nodd - 1
            case(LINE_2D_BOUNDARY_EAST)
                quarter_line(get_quarter_idx(nodd, 2), 1) = l - 1
                quarter_line(get_quarter_idx(nodd, 4), 1) = l - 1
                quarter_neighbour(get_quarter_idx(nodd, 2), 1) = nodu - 1
                quarter_neighbour(get_quarter_idx(nodd, 4), 1) = nodu - 1
            case(LINE_2D_BOUNDARY_SOUTH)
                quarter_line(get_quarter_idx(nodu, 1), 2) = l - 1
                quarter_line(get_quarter_idx(nodu, 2), 2) = l - 1
                quarter_neighbour(get_quarter_idx(nodu, 1), 2) = nodd - 1
                quarter_neighbour(get_quarter_idx(nodu, 2), 2) = nodd - 1
            case(LINE_2D_BOUNDARY_NORTH)
                quarter_line(get_quarter_idx(nodd, 3), 2) = l - 1
                quarter_line(get_quarter_idx(nodd, 4), 2) = l - 1
                quarter_neighbour(get_quarter_idx(nodd, 3), 2) = nodu - 1
                quarter_neighbour(get_quarter_idx(nodd, 4), 2) = nodu - 1
            end select
        enddo

    end subroutine set_quarter_admin


    function get_quarter_idx(nod, quadrant) result(idx)
        ! This function returns the quarter index based on the node and its quadrant index.
        !!!!!!!!!!!!!!!!
        !!  3  !!  4  !!
        !!!!!!!!!!!!!!!!
        !!  1  !!  2  !!
        !!!!!!!!!!!!!!!!
        integer, intent(in) :: nod
        integer, intent(in) :: quadrant
        integer :: idx

        idx = 4 * (nod - 1) + quadrant
    
    end function get_quarter_idx

end module m_cells
