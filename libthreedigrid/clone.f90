module m_clone

    use m_log
    
    use m_grid_utils, only : get_pix_corners
    
    implicit none

    contains

    subroutine find_active_clone_cells(n_cells, clone_array, cell_numbering, clone_numbering, clones_in_cell)
        ! ! Finding clone cells based on the clone array and renumbering the whole cells ! !
        
        use parameters, only : CLONE_NUMBERS
        integer, intent(inout) :: n_cells                           !! total number of active cells. in: without clone cells, out: with clone cells
        integer, intent(in) :: clone_array(:,:)                     !! Identifier of the clones within each host cell.
        integer, intent(inout) :: cell_numbering(:)                 !! New numbering for quadtree cells
        integer, intent(inout) :: clone_numbering(:)                !! New numbering for clone cells
        integer, intent(inout) :: clones_in_cell(:)                 !! Total number of clones in each host cell
        integer :: cell_counter, clone_counter, i, counter
        integer :: n_cells_new

        call print_unix('** INFO: Find active clone cells.')
        counter = 0 !! counter of the renumbering
        n_cells_new = n_cells
        do cell_counter = 1, n_cells !! Check which quadtree cells has clone cell(s) and how many
            if (clone_array(cell_counter, 1) >= 0) then !! If there are clone cells, count them in order
                do clone_counter = 1, CLONE_NUMBERS !! Count until the limit
                    if (clone_array(cell_counter, clone_counter) >= 0) then 
                        clone_numbering(clone_array(cell_counter, clone_counter) + 1) = counter !! to comply with python
                        counter = counter + 1
                        clones_in_cell(cell_counter) = clones_in_cell(cell_counter) + 1
                    end if                        
                end do
                if (clones_in_cell(cell_counter) == 1) then !! If there is one clone cell in a quadtree cell, discard it  (this will be moved to the filtering part in python)
                    clones_in_cell(cell_counter) = 0
                    clone_numbering(clone_array(cell_counter, 1) + 1) = -9999
                    cell_numbering(cell_counter) = counter - 1 !! (-1) to comply with python
                endif
            else !! If there are no clone cells, count the quadtree cells
                cell_numbering(cell_counter) = counter !! to comply with python
                counter = counter + 1
            end if
            if (clones_in_cell(cell_counter) > 0) then !! Ajdust the total number of the cells accordingly
                n_cells_new = n_cells_new + clones_in_cell(cell_counter) - 1
            end if
        end do
        n_cells = n_cells_new

        call print_unix('** INFO: Number of active quadtree and clone cells:', n_cells)

    end subroutine find_active_clone_cells

    subroutine make_clones(min_pix, area_mask, clone_mask, clones_in_cell, line, nodk, nodm, nodn, n_lines_u, n_lines_v, n_interclone_lines)
        !! Count the active (2D) lines created within the new cell system !!

        use parameters, only: U_DIR, V_DIR
        integer, intent(in) :: min_pix
        integer*1, intent(in) :: area_mask(:,:)
        integer, intent(in) :: clone_mask(:,:)                  !! Identifier of the clone cells according to the area mask
        integer, intent(in) :: clones_in_cell(:)                !! Updated identifier of the clones within each host cell (does not match the clone_mask)
        integer, intent(in) :: line(:,:)                        !! Old Line administration
        integer, intent(in) :: nodk(:)                          !! Old refinement level administration
        integer, intent(in) :: nodm(:)                          !! Old row number administration
        integer, intent(in) :: nodn(:)                          !! Old column number administration
        integer, intent(inout) :: n_lines_u                     !! Total number of 2D line in u direction. in: old number, out: new number including clones
        integer, intent(inout) :: n_lines_v                     !! Total number of 2D line in v direction. in: old number, out: new number including clones
        integer, intent(inout) :: n_interclone_lines            !! Total number of line linking 2 clones within the same quadtree cell.

        integer :: start_u, end_u, start_v, end_v

        start_u = 1
        end_u = n_lines_u
        start_v = n_lines_u + 1
        end_v = n_lines_u + n_lines_v

        call count_clone_lines(U_DIR, start_u, end_u, n_lines_u, min_pix, area_mask, clone_mask, clones_in_cell, line, nodk, nodm, nodn)
        
        call count_clone_lines(V_DIR, start_v, end_v, n_lines_v, min_pix, area_mask, clone_mask, clones_in_cell, line, nodk, nodm, nodn)

        call count_interclone_lines(n_interclone_lines, clones_in_cell)

        call print_unix('** INFO: Number of 2D lines in u- and v-dir:', n_lines_u, n_lines_v)
        call print_unix('** INFO: Number of inter-clone lines:', n_interclone_lines)

    end subroutine

    subroutine count_clone_lines(direction, start_l, end_l, n_line, min_pix, area_mask, clone_mask, clones_in_cell, line, nodk, nodm, nodn)

        use parameters, only: U_DIR, V_DIR
        integer, intent(in) :: direction
        integer, intent(in) :: start_l
        integer, intent(in) :: end_l
        integer, intent(inout) :: n_line
        integer, intent(in) :: min_pix
        integer*1, intent(in) :: area_mask(:,:)
        integer, intent(in) :: clone_mask(:,:)                !! Identifier of the clone cells according to the area mask
        integer, intent(in) :: clones_in_cell(:)              !! total number of clones in each cell
        integer, intent(in) :: line(:,:)
        integer, intent(in) :: nodk(:)
        integer, intent(in) :: nodm(:)
        integer, intent(in) :: nodn(:)

        integer :: host_1, host_2, k_host_1, k_host_2, clone_1, clone_2
        integer :: pixel_i, pixel_j, i0, i1, j0, j1, nl
        integer :: num_pix, n_line_new, pixel_no

        n_line_new = n_line 
        do nl = start_l, end_l
            host_1 = line(nl,1) + 1
            host_2 = line(nl,2) + 1
            if (clones_in_cell(host_1) > 0 .or. clones_in_cell(host_2) > 0) then !! It means that the quadtree flowline must be deleted
                k_host_1 = nodk(host_1)
                k_host_2 = nodk(host_2)

                num_pix = min_pix * 2**(min(k_host_1,k_host_2)-1) !! the minimum number of the pixels within two adjacent cells

                n_line_new = n_line_new - 1   !! One quadtree flowline is deleted
                call get_pix_corners(k_host_1, nodm(host_1), nodn(host_1), min_pix, i0, i1, j0, j1)
                !! Pick the (i,j) of the pixel at the bottom-right corner of the left cell
                if (direction == U_DIR) then
                    pixel_i = i1
                    pixel_j = j0
                    call get_pix_corners(k_host_2, nodm(host_2), nodn(host_2), min_pix, i0, i1, j0, j1)
                    pixel_j = max(j0, pixel_j)   !! Helps with the correct refinement level
                    pixel_i = min(i1, pixel_i)   !! Helps with the correct selection of the left cell
                    clone_1 = clone_mask(pixel_i, pixel_j)
                    clone_2 = clone_mask(pixel_i + 1, pixel_j)
                elseif (direction == V_DIR) then
                    pixel_i = i0
                    pixel_j = j1
                    call get_pix_corners(k_host_2, nodm(host_2), nodn(host_2), min_pix, i0, i1, j0, j1)
                    pixel_i = max(i0, pixel_i)   !! Helps with the correct refinement level
                    pixel_j = min(j1, pixel_j)   !! Helps with the correct selection of the left cell
                    clone_1 = clone_mask(pixel_i, pixel_j)
                    clone_2 = clone_mask(pixel_i, pixel_j + 1)
                endif
                !! Now check whether or not there are at least one pair of active pixels btw 2 cells
                pixel_no = 0
                call find_active_lines(direction, pixel_i, pixel_j, min_pix, clone_1, clone_2, num_pix, area_mask, clone_mask, n_line_new, pixel_no)   !! Only counting the lines         
            end if
        end do
        
        n_line = n_line_new

    end subroutine count_clone_lines

    subroutine count_interclone_lines(n_lines, clones_in_cell)
        
        integer, intent(inout) :: n_lines  !! number of interclone lines
        integer, intent(in) :: clones_in_cell(:)

        integer :: n_cell

        do n_cell = 1, size(clones_in_cell)
            if (clones_in_cell(n_cell) == 2) then
                n_lines = n_lines + 1
            ! elseif (clones_in_cell(n_cell) == 2) then !!!! TODO: We have to cover upto 4 clone cells in a quadtree cell
            endif
        enddo

    end subroutine count_interclone_lines

    subroutine set_lines_with_clones(n_lines_u, n_lines_v, n_interclone_lines, min_pix, area_mask, clone_mask, clone_array, clones_in_cell, cell_numbering, clone_numbering, line, nodk, nodm, nodn, line_new, cross_pix_new, cross_pix)

        use parameters, only: U_DIR, V_DIR
        integer, intent(inout) :: n_lines_u
        integer, intent(inout) :: n_lines_v
        integer, intent(in) :: n_interclone_lines
        integer, intent(in) :: min_pix
        integer*1, intent(in) :: area_mask(:,:)
        integer, intent(in) :: clone_mask(:,:)             !! Identifier of the clone cells according to the area mask
        integer, intent(in) :: clone_array(:,:)             !! Identifier of the clone cells according to the area mask
        integer, intent(in) :: clones_in_cell(:)              !! Updated identifier of the clones within each host cell (does not match the clone_mask)
        integer, intent(in) :: cell_numbering(:)             !! New numbering list for quadtree cells
        integer, intent(in) :: clone_numbering(:)            !! New numbering list for clone cells
        integer, intent(in) :: line(:,:)
        integer, intent(in) :: nodk(:)
        integer, intent(in) :: nodm(:)
        integer, intent(in) :: nodn(:)
        integer, intent(inout) :: line_new(:,:)
        integer, intent(inout) :: cross_pix_new(:,:)
        integer, intent(in) :: cross_pix(:,:)

        integer :: start_u, end_u, start_v, end_v, counter

        call print_unix('** INFO: Resetting line attribute...')

        start_u = 1
        end_u = n_lines_u
        start_v = n_lines_u + 1
        end_v = n_lines_u + n_lines_v

        counter = 0

        call set_2D_lines_with_clones(U_DIR, start_u, end_u, n_lines_u, min_pix, area_mask, clone_mask, clones_in_cell, cell_numbering, clone_numbering, line, nodk, nodm, nodn, counter, line_new, cross_pix_new, cross_pix)
        
        call set_2D_lines_with_clones(V_DIR, start_v, end_v, n_lines_v, min_pix, area_mask, clone_mask, clones_in_cell, cell_numbering, clone_numbering, line, nodk, nodm, nodn, counter, line_new, cross_pix_new, cross_pix)
        
        call set_interclone_lines(clones_in_cell, counter, line_new, clone_numbering, clone_array)

    end subroutine

    subroutine set_2d_lines_with_clones(direction, start_l, end_l, n_line_new, min_pix, area_mask, clone_mask, clones_in_cell, cell_numbering, clone_numbering, line, nodk, nodm, nodn, l_counter, line_new, cross_pix_new, cross_pix)

        use parameters, only: U_DIR, V_DIR
        integer, intent(in) :: direction
        integer, intent(in) :: start_l
        integer, intent(in) :: end_l
        integer, intent(inout) :: n_line_new
        integer, intent(in) :: min_pix
        integer*1, intent(in) :: area_mask(:,:)
        integer, intent(in) :: clone_mask(:,:)             !! Identifier of the clone cells according to the area mask
        integer, intent(in) :: clones_in_cell(:)              !! Updated identifier of the clones within each host cell (does not match the clone_mask)
        integer, intent(in) :: cell_numbering(:)             !! New numbering list for quadtree cells
        integer, intent(in) :: clone_numbering(:)            !! New numbering list for clone cells
        integer, intent(in) :: line(:,:)
        integer, intent(in) :: nodk(:)
        integer, intent(in) :: nodm(:)
        integer, intent(in) :: nodn(:)
        integer, intent(inout) :: l_counter
        integer, intent(inout) :: line_new(:,:)
        integer, intent(inout) :: cross_pix_new(:,:)
        integer, intent(in) :: cross_pix(:,:)

        integer :: host_1, host_2, k_host_1, k_host_2, m_host_1, n_host_1, clone_1, clone_2
        integer :: pixel_i, pixel_j, i0, i1, j0, j1, nl, pixel_no, num_pix
    
        do nl = start_l, end_l
            host_1 = line(nl,1) + 1
            host_2 = line(nl,2) + 1
            if (clones_in_cell(host_1) > 0 .or. clones_in_cell(host_2) > 0) then
                k_host_1 = nodk(host_1)
                k_host_2 = nodk(host_2)
                num_pix = min_pix * 2**(min(k_host_1,k_host_2)-1) !! temporary

                call get_pix_corners(k_host_1, nodm(host_1), nodn(host_1), min_pix, i0, i1, j0, j1)
                if (direction == U_DIR) then
                    pixel_i = i1
                    pixel_j = j0
                    call get_pix_corners(k_host_2, nodm(host_2), nodn(host_2), min_pix, i0, i1, j0, j1)
                    pixel_j = max(j0, pixel_j)
                    pixel_i = min(i1, pixel_i)
                    clone_1 = clone_mask(pixel_i, pixel_j)
                    clone_2 = clone_mask(pixel_i + 1, pixel_j)
                elseif (direction == V_DIR) then
                    pixel_i = i0
                    pixel_j = j1
                    call get_pix_corners(k_host_2, nodm(host_2), nodn(host_2), min_pix, i0, i1, j0, j1)
                    pixel_i = max(i0, pixel_i)
                    pixel_j = min(j1, pixel_j)
                    clone_1 = clone_mask(pixel_i, pixel_j)
                    clone_2 = clone_mask(pixel_i, pixel_j + 1)
                endif

                pixel_no = 0
                call find_active_lines(direction, pixel_i, pixel_j, min_pix, clone_1, clone_2, num_pix, area_mask, clone_mask, n_line_new, pixel_no, &
                                            line_new, host_1, host_2, l_counter, clone_numbering, cell_numbering, cross_pix_new)   !! Now rewiring the line administration
            else
                l_counter = l_counter + 1
                line_new(l_counter, 1) = cell_numbering(line(nl,1) + 1)
                line_new(l_counter, 2) = cell_numbering(line(nl,2) + 1)
                cross_pix_new(l_counter,:) = cross_pix(nl,:)
            end if
        end do
    end subroutine set_2d_lines_with_clones
    
    recursive subroutine find_active_lines(direction, pixel_i, pixel_j, min_pix, clone_1, clone_2, num_pix, area_mask, clone_mask, tot_number_lines, pixel_no, &
                                                line_new, host_1, host_2, l_counter, clone_numbering, cell_numbering, cross_pix_new)
        
        use parameters, only: U_DIR, V_DIR
        integer, intent(in) :: direction
        integer, intent(in) :: pixel_i
        integer, intent(in) :: pixel_j
        integer, intent(in) :: min_pix
        integer, intent(in) :: clone_1
        integer, intent(in) :: clone_2
        integer, intent(in) :: num_pix
        integer*1, intent(in) :: area_mask(:,:)
        integer, intent(in) :: clone_mask(:,:)
        integer, intent(inout) :: tot_number_lines
        integer, intent(inout) :: pixel_no
        integer, intent(inout), optional :: line_new(:,:)
        integer, intent(in), optional :: host_1
        integer, intent(in), optional :: host_2
        integer, intent(inout), optional :: l_counter
        integer, intent(in), optional :: clone_numbering(:)
        integer, intent(in), optional :: cell_numbering(:)
        integer, intent(inout), optional :: cross_pix_new(:,:)

        integer :: counter, pixel, new_pixel
        integer :: new_clone_1, new_clone_2

        if (direction == U_DIR) then  !! u-direction
            pixel = pixel_j + pixel_no
            do counter = pixel, pixel_j + num_pix - 1  !! find the common border
                if (clone_mask(pixel_i, counter) == clone_1 .and. clone_mask(pixel_i+1, counter) == clone_2) then
                    pixel_no = pixel_no + 1
                else
                    exit
                endif
            enddo
            new_pixel = pixel_j + pixel_no
            if ((new_pixel - pixel) > 1) then   !! 1 is the minimum number of pixels through which a flowline can be created
                if (any(minval(area_mask(pixel_i:pixel_i+1,pixel:new_pixel-1), 1) > 0)) then  !! if there is a minimum of a pair of pixel with data, make the flowline
                    if (present(line_new)) then
                        l_counter = l_counter + 1
                        cross_pix_new(l_counter,:) = (/pixel_i, pixel-1, pixel_i, new_pixel-1/)
                        if (clone_1 >= 0 .and. clone_2 >= 0) then
                            line_new(l_counter, 1) = clone_numbering(clone_1+1)
                            line_new(l_counter, 2) = clone_numbering(clone_2+1)
                        elseif (clone_1 >= 0 .and. clone_2 < 0) then
                            line_new(l_counter, 1) = clone_numbering(clone_1+1)
                            line_new(l_counter, 2) = cell_numbering(host_2)
                        elseif (clone_1 < 0 .and. clone_2 >= 0 ) then
                            line_new(l_counter, 1) = cell_numbering(host_1)
                            line_new(l_counter, 2) = clone_numbering(clone_2+1)
                        endif
                    else
                        tot_number_lines = tot_number_lines + 1
                    endif
                endif
            endif
            if (pixel_no < num_pix) then
                do counter = new_pixel, pixel_j + num_pix - 1  !! Check if both clone_1&2 change
                    new_clone_1 = clone_mask(pixel_i, counter)
                    new_clone_2 = clone_mask(pixel_i + 1, counter)
                    if (clone_1 < 0 .and. new_clone_2 /= clone_2) then
                        exit
                    elseif (new_clone_1 /= clone_1 .and. clone_2 < 0) then
                        exit
                    elseif (new_clone_1 /= clone_1 .and. new_clone_2 /= clone_2) then
                        exit
                    else
                        pixel_no = pixel_no + 1
                    endif
                enddo

                ! new_clone_1 = clone_mask(pixel_i, new_pixel)
                ! new_clone_2 = clone_mask(pixel_i + 1, new_pixel)
                
                if (present(line_new)) then
                    call find_active_lines(direction, pixel_i, pixel_j, min_pix, new_clone_1, new_clone_2, num_pix, area_mask, clone_mask, tot_number_lines, pixel_no, &
                                                line_new, host_1, host_2, l_counter, clone_numbering, cell_numbering, cross_pix_new)
                else
                    call find_active_lines(direction, pixel_i, pixel_j, min_pix, new_clone_1, new_clone_2, num_pix, area_mask, clone_mask, tot_number_lines, pixel_no)
                endif
            endif

        elseif (direction == V_DIR) then  !! v-direction
            pixel = pixel_i + pixel_no
            do counter = pixel, pixel_i + num_pix - 1  !! find the common border
                if (clone_mask(counter, pixel_j) == clone_1 .and. clone_mask(counter, pixel_j+1) == clone_2) then
                    pixel_no = pixel_no + 1
                else
                    exit
                endif
            enddo
            new_pixel = pixel_i + pixel_no
            if ((new_pixel - pixel) > 1) then
                if (any(minval(area_mask(pixel:new_pixel-1,pixel_j:pixel_j+1), 2) > 0)) then
                    if (present(line_new)) then
                        l_counter = l_counter + 1
                        cross_pix_new(l_counter,:) = (/pixel-1, pixel_j, new_pixel-1, pixel_j/)
                        if (clone_1 >= 0 .and. clone_2 >= 0) then
                            line_new(l_counter, 1) = clone_numbering(clone_1+1)
                            line_new(l_counter, 2) = clone_numbering(clone_2+1)
                        elseif (clone_1 >= 0 .and. clone_2 < 0) then
                            line_new(l_counter, 1) = clone_numbering(clone_1+1)
                            line_new(l_counter, 2) = cell_numbering(host_2)
                        elseif (clone_1 < 0 .and. clone_2 >= 0) then
                            line_new(l_counter, 1) = cell_numbering(host_1)
                            line_new(l_counter, 2) = clone_numbering(clone_2+1)
                        endif
                    else
                        tot_number_lines = tot_number_lines + 1
                    endif
                endif
            endif
            if (pixel_no < num_pix) then
                do counter = new_pixel, pixel_i + num_pix - 1 !! Check if both clone_1&2 change
                    new_clone_1 = clone_mask(counter, pixel_j)
                    new_clone_2 = clone_mask(counter, pixel_j + 1)
                    if (clone_1 < 0 .and. new_clone_2 /= clone_2) then
                        exit
                    elseif (new_clone_1 /= clone_1 .and. clone_2 < 0) then
                        exit
                    elseif (new_clone_1 /= clone_1 .and. new_clone_2 /= clone_2) then
                        exit
                    else
                        pixel_no = pixel_no + 1
                    endif
                enddo

                ! new_clone_1 = clone_mask(new_pixel, pixel_j)
                ! new_clone_2 = clone_mask(new_pixel, pixel_j + 1)
                
                if (present(line_new)) then
                    call find_active_lines(direction, pixel_i, pixel_j, min_pix, new_clone_1, new_clone_2, num_pix, area_mask, clone_mask, tot_number_lines, pixel_no, &
                                                line_new, host_1, host_2, l_counter, clone_numbering, cell_numbering, cross_pix_new)
                else
                    call find_active_lines(direction, pixel_i, pixel_j, min_pix, new_clone_1, new_clone_2, num_pix, area_mask, clone_mask, tot_number_lines, pixel_no)
                endif
            endif

        endif

    end subroutine find_active_lines

    subroutine set_interclone_lines(clones_in_cell, l_counter, line_new, clone_numbering, clone_array)
        
        integer, intent(in) :: clones_in_cell(:)
        integer, intent(inout):: l_counter
        integer, intent(inout) :: line_new(:,:)
        integer, intent(in) :: clone_numbering(:)
        integer, intent(in) :: clone_array(:,:)

        integer :: n_cell

        do n_cell = 1, size(clones_in_cell)
            if (clones_in_cell(n_cell) == 2) then
                l_counter = l_counter + 1
                line_new(l_counter, 1) = clone_numbering(clone_array(n_cell, 1) + 1)
                line_new(l_counter, 2) = clone_numbering(clone_array(n_cell, 2) + 1)
                ! elseif (clones_in_cell(n_cell) == 2) then  !!! TODO: We will cover up to 4 clones per quadtree cell
            endif
        enddo

    end subroutine set_interclone_lines

    subroutine reset_nod_parameters(clone_array, nodk, nodm, nodn, nodk_new, nodm_new, nodn_new,&
                                    bounds, coords, pixel_coords, bounds_new, coords_new, pixel_coords_new)

        use parameters, only : CLONE_NUMBERS
        integer, intent(in) :: clone_array(:,:)
        integer, intent(in) :: nodk(:)
        integer, intent(in) :: nodm(:)
        integer, intent(in) :: nodn(:)
        integer, intent(inout) :: nodk_new(:)
        integer, intent(inout) :: nodm_new(:)
        integer, intent(inout) :: nodn_new(:)
        double precision, intent(in) :: bounds(:,:)
        double precision, intent(in) :: coords(:,:)
        integer, intent(in) :: pixel_coords(:,:)
        double precision, intent(inout) :: bounds_new(:,:)
        double precision, intent(inout) :: coords_new(:,:)
        integer, intent(inout) :: pixel_coords_new(:,:)

        integer :: cell_counter, clone_counter, counter

        call print_unix('** INFO: Reset nodes attributes...')

        counter = 0
        do cell_counter = 1, size(nodk)
            if (clone_array(cell_counter, 1) >= 0) then
                do clone_counter = 1, CLONE_NUMBERS !! Count until the limit
                    if (clone_array(cell_counter, clone_counter) >= 0) then
                        counter = counter + 1
                        nodk_new(counter) = nodk(cell_counter)
                        nodm_new(counter) = nodm(cell_counter)
                        nodn_new(counter) = nodn(cell_counter)
                        bounds_new(counter, :) = bounds(cell_counter, :)
                        coords_new(counter, :) = coords(cell_counter, :)
                        pixel_coords_new(counter, :) = pixel_coords(cell_counter, :)
                    endif
                enddo
            else
                counter = counter + 1
                nodk_new(counter) = nodk(cell_counter)
                nodm_new(counter) = nodm(cell_counter)
                nodn_new(counter) = nodn(cell_counter)
                bounds_new(counter, :) = bounds(cell_counter, :)
                coords_new(counter, :) = coords(cell_counter, :)
                pixel_coords_new(counter, :) = pixel_coords(cell_counter, :)
            endif
        enddo
        
    end subroutine reset_nod_parameters

    subroutine set_line_coords_new(line_number, line, nodk, nodm, nodn, node_coord, node_bound, min_pix, clone_mask, clone_numbering, clone_centroid, clone_polygon, line_coords, centroids, polygons)

        integer, intent(in) :: line_number
        integer, intent(in) :: line(:,:)
        integer, intent(in) :: nodk(:)
        integer, intent(in) :: nodm(:)
        integer, intent(in) :: nodn(:)
        double precision, intent(in) :: node_coord(:,:)
        double precision, intent(in) :: node_bound(:,:,:)
        integer, intent(in) :: min_pix
        integer, intent(in) :: clone_mask(:,:)
        integer, intent(in) :: clone_numbering(:)
        double precision, intent(in) :: clone_centroid(:,:)
        double precision, intent(in) :: clone_polygon(:,:,:)
        double precision, intent(inout) :: line_coords(:,:)      ! center coords of comp cell at two ends
        double precision, intent(inout) :: centroids(:,:)        ! centroids of all the nodes (with clone cells)
        double precision, intent(inout) :: polygons(:,:,:)       ! bounds of all the cell polygons (with clone cells)
        integer :: line_counter, cell_no, find_index
        integer :: i0, j0, i1, j1
        integer*1, allocatable :: centroid_key(:)

        allocate(centroid_key(size(centroids, 1)))
        centroid_key = 0

        do line_counter = 1, line_number
            cell_no = line(line_counter, 1) + 1
            call get_pix_corners(nodk(cell_no), nodm(cell_no), nodn(cell_no), min_pix, i0, i1, j0, j1)
            if (clone_mask(i0, j0) >= 0) then
                do find_index = 1, size(clone_numbering,1)
                    if (clone_numbering(find_index) == cell_no - 1) then
                        exit
                    endif
                enddo
                line_coords(line_counter, 1) = clone_centroid(find_index,1)
                line_coords(line_counter, 2) = clone_centroid(find_index,2)
                if (centroid_key(cell_no) == 0) then
                    centroids(cell_no, :) = clone_centroid(find_index, :)
                    polygons(cell_no, :, :) = clone_polygon(find_index, :, :)
                    centroid_key(cell_no) = 1
                endif
            else
                line_coords(line_counter, 1) = node_coord(cell_no,1)
                line_coords(line_counter, 2) = node_coord(cell_no,2)
                if (centroid_key(cell_no) == 0) then
                    centroids(cell_no, :) = node_coord(cell_no, :)
                    polygons(cell_no, :, :) = node_bound(cell_no, :, :)
                    centroid_key(cell_no) = 1
                endif
            endif
            cell_no = line(line_counter, 2) + 1
            call get_pix_corners(nodk(cell_no), nodm(cell_no), nodn(cell_no), min_pix, i0, i1, j0, j1)
            if (clone_mask(i0, j0) >= 0) then
                do find_index = 1, size(clone_numbering,1)
                    if (clone_numbering(find_index) == cell_no - 1) then
                        exit
                    endif
                enddo
                line_coords(line_counter, 3) = clone_centroid(find_index,1)
                line_coords(line_counter, 4) = clone_centroid(find_index,2)
                if (centroid_key(cell_no) == 0) then
                    centroids(cell_no, :) = clone_centroid(find_index, :)
                    polygons(cell_no, :, :) = clone_polygon(find_index, :, :)
                    centroid_key(cell_no) = 1
                endif
            else
                line_coords(line_counter, 3) = node_coord(cell_no,1)
                line_coords(line_counter, 4) = node_coord(cell_no,2)
                if (centroid_key(cell_no) == 0) then
                    centroids(cell_no, :) = node_coord(cell_no, :)
                    polygons(cell_no, :, :) = node_bound(cell_no, :, :)
                    centroid_key(cell_no) = 1
                endif
            endif
        enddo

    end subroutine set_line_coords_new

    subroutine set_quad_idx(quad_idx, nodk, nodm, nodn, n_cells)

        use m_grid_utils, only : get_lg_corners
        
        integer, intent(inout) :: quad_idx(:,:)
        integer, intent(in) :: nodk(:)
        integer, intent(in) :: nodm(:)
        integer, intent(in) :: nodn(:)
        integer, intent(in) ::n_cells

        integer :: quad_idx_new(size(quad_idx, 1), size(quad_idx, 2)), mn(4)
        integer :: cell, k, m, n

        quad_idx_new = 0.0d0

        do cell = 1, n_cells
            k = nodk(cell)
            m = nodm(cell)
            n = nodn(cell)
            mn = get_lg_corners(k, m, n)
            if (all(quad_idx_new(mn(1):mn(3), mn(2):mn(4)) == 0.0d0)) then
                quad_idx_new(mn(1):mn(3), mn(2):mn(4))= cell
            endif
        enddo

        quad_idx = quad_idx_new

    end subroutine set_quad_idx

    subroutine set_visualization_cell_bounds(bounds, cell_bounds)

        !! Propogate the coords of bottom-left and top-right corners to all corners
        double precision, intent(in) :: bounds(:,:)                   ! corner coords of quadtree cells at bottom-left and top-right
        double precision, intent(inout) :: cell_bounds(:,:,:)         ! all corner coords of all the quadtree cells

        integer :: cell

        do cell = 1, size(bounds, 1)
            cell_bounds(cell, 1, :) = bounds(cell, 1:2)
            cell_bounds(cell, 2, 1) = bounds(cell, 1)
            cell_bounds(cell, 2, 2) = bounds(cell, 4)
            cell_bounds(cell, 3, :) = bounds(cell, 3:4)
            cell_bounds(cell, 4, 1) = bounds(cell, 3)
            cell_bounds(cell, 4, 2) = bounds(cell, 2)
            cell_bounds(cell, 5, :) = bounds(cell, 1:2)
        enddo

    end subroutine set_visualization_cell_bounds


end module