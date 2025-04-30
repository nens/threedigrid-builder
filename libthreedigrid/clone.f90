module m_clone

    use m_grid_utils, only : get_pix_corners
    
    implicit none

    contains

    subroutine find_active_clone_cells(n_cells, clone_array, cell_numbering, clone_numbering, clones_in_cell)  !!, n_clone_cells, clones_in_cell
        ! ! Finding clone cells based on the area mask specifically defined for clone cells
        
        use m_grid_utils, only : get_pix_corners
        use parameters, only : CLONE_NUMBERS
        integer, intent(inout) :: n_cells                  !! total number of active cells. in: without clone cells, out: with clone cells
        integer, intent(inout) :: clone_array(:,:)         !! Identifier of the clones within each host cell. in: matches the clone_mask, out: matches the new numbering
        integer, intent(inout) :: cell_numbering(:)        !! New numbering for quadtree cells
        integer, intent(inout) :: clone_numbering(:)       !! New numbering for clone cells
        integer, intent(inout) :: clones_in_cell(:)                 !! Total number of clones in each host cell
        ! integer :: clones_in_cell(n_cells)                 !! Total number of clones in each host cell
        integer :: cell_counter, clone_counter, i, counter
        integer :: n_cells_new

        counter = 0 !! counter of the renumbering
        n_cells_new = n_cells
        do cell_counter = 1, n_cells !! Check which quadtree cells has clone cell(s) and how many
            if (clone_array(cell_counter, 1) >= 0) then !! If there are clone cells, count them in order
                do clone_counter = 1, CLONE_NUMBERS !! Count until the limit
                    if (clone_array(cell_counter, clone_counter) >= 0) then 
                        counter = counter + 1
                        clone_numbering(clone_array(cell_counter, clone_counter) + 1) = counter
                        clones_in_cell(cell_counter) = clones_in_cell(cell_counter) + 1
                    end if                        
                end do
                if (clones_in_cell(cell_counter) == 1) then !! If there is one clone cell in a quadtree cell, discard it
                    clones_in_cell(cell_counter) = 0
                    clone_numbering(clone_array(cell_counter, 1) + 1) = 0
                    cell_numbering(cell_counter) = counter
                endif
            else !! If there are no clone cells, count the quadtree cells
                counter = counter + 1
                cell_numbering(cell_counter) = counter
            end if
            if (clones_in_cell(cell_counter) > 0) then !! Ajdust the total number of the cells accordingly
                n_cells_new = n_cells_new + clones_in_cell(cell_counter) - 1
            end if
        end do
        n_cells = n_cells_new

    end subroutine find_active_clone_cells

    subroutine count_clone_lines(direction, start_l, end_l, n_line, lgrmin, area_mask, clone_mask, clones_in_cell, cell_numbering, clone_numbering, line, nodk, nodm, nodn)

        use parameters, only: U_DIR, V_DIR
        integer, intent(in) :: direction
        integer, intent(in) :: start_l
        integer, intent(in) :: end_l
        integer, intent(inout) :: n_line
        integer, intent(in) :: lgrmin
        integer*1, intent(in) :: area_mask(:,:)
        integer, intent(in) :: clone_mask(:,:)             !! Identifier of the clone cells according to the area mask
        integer, intent(in) :: clones_in_cell(:)              !! Updated identifier of the clones within each host cell (does not match the clone_mask)
        integer, intent(in) :: cell_numbering(:)             !! New numbering list for quadtree cells
        integer, intent(in) :: clone_numbering(:)            !! New numbering list for clone cells
        integer, intent(inout) :: line(:,:)
        integer, intent(in) :: nodk(:)
        integer, intent(in) :: nodm(:)
        integer, intent(in) :: nodn(:)

        integer :: host_1, host_2, k_host_1, k_host_2, m_host_1, n_host_1, clone_1, clone_2
        integer :: pixel_i, pixel_j, i0, i1, j0, j1, l_counter, nl
        integer :: num_pix    !! temp
        integer :: n_line_new    !! temp
        integer :: pixel_no
        integer, allocatable :: cell_order(:)      ! administration of column/row no of cell wrt the direction
        

        n_line_new = n_line 
        do nl = start_l, end_l   !!size(line_new) n_line_u, n_line_v
            host_1 = line(nl,1) + 1
            host_2 = line(nl,2) + 1
            if (clones_in_cell(host_1) > 0 .or. clones_in_cell(host_2) > 0) then !! So it means that the quadtree flowline must be deleted
                k_host_1 = nodk(host_1)
                k_host_2 = nodk(host_2)

                num_pix = lgrmin * 2**(min(k_host_1,k_host_2)-1) !! temporary

                n_line_new = n_line_new - 1   !! if no refinement or refinement at the left-side, there is one flowline
                call get_pix_corners(k_host_1, nodm(host_1), nodn(host_1), lgrmin, i0, i1, j0, j1)
                if (direction == U_DIR) then
                    pixel_i = i1
                    pixel_j = j0
                    call get_pix_corners(k_host_2, nodm(host_2), nodn(host_2), lgrmin, i0, i1, j0, j1)
                    pixel_j = max(j0, pixel_j)
                    pixel_i = min(i1, pixel_i)
                    clone_1 = clone_mask(pixel_i, pixel_j)
                    clone_2 = clone_mask(pixel_i + 1, pixel_j)
                elseif (direction == V_DIR) then
                    pixel_i = i0
                    pixel_j = j1
                    call get_pix_corners(k_host_2, nodm(host_2), nodn(host_2), lgrmin, i0, i1, j0, j1)
                    pixel_i = max(i0, pixel_i)
                    pixel_j = min(j1, pixel_j)
                    clone_1 = clone_mask(pixel_i, pixel_j)
                    clone_2 = clone_mask(pixel_i, pixel_j + 1)
                endif
                pixel_no = 0
                call find_active_clone_lines(direction, pixel_i, pixel_j, lgrmin, clone_1, clone_2, num_pix, area_mask, clone_mask, n_line_new, pixel_no)   !! Only counting the lines         
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
            ! elseif (clones_in_cell(n_cell) == 2) then
            endif
        enddo

    end subroutine count_interclone_lines

    subroutine find_clone_lines(direction, start_l, end_l, n_line_new, lgrmin, area_mask, clone_mask, clones_in_cell, cell_numbering, clone_numbering, line, nodk, nodm, nodn, l_counter, line_new)

        use parameters, only: U_DIR, V_DIR
        integer, intent(in) :: direction
        integer, intent(in) :: start_l
        integer, intent(in) :: end_l
        integer, intent(inout) :: n_line_new
        integer, intent(in) :: lgrmin
        integer*1, intent(in) :: area_mask(:,:)
        integer, intent(in) :: clone_mask(:,:)             !! Identifier of the clone cells according to the area mask
        integer, intent(in) :: clones_in_cell(:)              !! Updated identifier of the clones within each host cell (does not match the clone_mask)
        integer, intent(in) :: cell_numbering(:)             !! New numbering list for quadtree cells
        integer, intent(in) :: clone_numbering(:)            !! New numbering list for clone cells
        integer, intent(inout) :: line(:,:)
        integer, intent(in) :: nodk(:)
        integer, intent(in) :: nodm(:)
        integer, intent(in) :: nodn(:)
        integer, intent(inout) :: l_counter
        integer, intent(inout) :: line_new(:,:)

        integer :: host_1, host_2, k_host_1, k_host_2, m_host_1, n_host_1, clone_1, clone_2
        integer :: pixel_i, pixel_j, i0, i1, j0, j1, nl, pixel_no, num_pix
        integer, allocatable :: cell_order(:)      ! administration of column/row no of cell wrt the direction
    
        do nl = start_l, end_l
            host_1 = line(nl,1) + 1
            host_2 = line(nl,2) + 1
            if (clones_in_cell(host_1) > 0 .or. clones_in_cell(host_2) > 0) then
                k_host_1 = nodk(host_1)
                k_host_2 = nodk(host_2)
                num_pix = lgrmin * 2**(min(k_host_1,k_host_2)-1) !! temporary

                call get_pix_corners(k_host_1, nodm(host_1), nodn(host_1), lgrmin, i0, i1, j0, j1)
                if (direction == U_DIR) then
                    pixel_i = i1
                    pixel_j = j0
                    call get_pix_corners(k_host_2, nodm(host_2), nodn(host_2), lgrmin, i0, i1, j0, j1)
                    pixel_j = max(j0, pixel_j)
                    pixel_i = min(i1, pixel_i)
                    clone_1 = clone_mask(pixel_i, pixel_j)
                    clone_2 = clone_mask(pixel_i + 1, pixel_j)
                elseif (direction == V_DIR) then
                    pixel_i = i0
                    pixel_j = j1
                    call get_pix_corners(k_host_2, nodm(host_2), nodn(host_2), lgrmin, i0, i1, j0, j1)
                    pixel_i = max(i0, pixel_i)
                    pixel_j = min(j1, pixel_j)
                    clone_1 = clone_mask(pixel_i, pixel_j)
                    clone_2 = clone_mask(pixel_i, pixel_j + 1)
                endif

                pixel_no = 0
                call find_active_clone_lines(direction, pixel_i, pixel_j, lgrmin, clone_1, clone_2, num_pix, area_mask, clone_mask, n_line_new, pixel_no, &
                                            line_new, host_1, host_2, l_counter, clone_numbering, cell_numbering)   !! Now rewiring the line administration
            else
                l_counter = l_counter + 1
                line_new(l_counter, 1) = cell_numbering(line(nl,1) + 1)
                line_new(l_counter, 2) = cell_numbering(line(nl,2) + 1)
            end if
        end do
    end subroutine find_clone_lines
    
    recursive subroutine find_active_clone_lines(direction, pixel_i, pixel_j, min_pix, clone_1, clone_2, num_pix, area_mask, clone_mask, tot_number_lines, pixel_no, &
                                                line_new, host_1, host_2, l_counter, clone_numbering, cell_numbering)
        
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

        integer :: counter, pixel, new_pixel
        integer :: new_clone_1, new_clone_2
        integer :: new_pixel1, new_pixel2 !!temp

        integer :: temp

        temp = 0
        if (direction == U_DIR) then  !! u-direction
            pixel = pixel_j + pixel_no
            do counter = pixel, pixel + num_pix - 1  !! find the common border
                if (clone_mask(pixel_i, counter) == clone_1 .and. clone_mask(pixel_i+1, counter) == clone_2) then
                    pixel_no = pixel_no + 1
                else
                    exit
                endif
            enddo
            new_pixel = pixel_j + pixel_no !! I can also use counter instead of new_pixel
            if ((new_pixel - pixel) > 1) then
                if (any(minval(area_mask(pixel_i:pixel_i+1,pixel:new_pixel-1), 1) > 0)) then  !! if there is a minimum of a pair of pixel with data, make the flowline
                    if (present(line_new)) then
                        l_counter = l_counter + 1
                        if (clone_1 >= 0 .and. clone_2 >= 0) then
                            if (clone_numbering(clone_1 + 1) > 0 .and. clone_numbering(clone_2 + 1) > 0) then
                                line_new(l_counter, 1) = clone_numbering(clone_1+1)
                                line_new(l_counter, 2) = clone_numbering(clone_2+1)
                            elseif (clone_numbering(clone_1 + 1) == 0 .and. clone_numbering(clone_2 + 1) > 0) then
                                line_new(l_counter, 1) = cell_numbering(host_1)
                                line_new(l_counter, 2) = clone_numbering(clone_2+1)
                            elseif (clone_numbering(clone_1 + 1) > 0 .and. clone_numbering(clone_2 + 1) == 0) then
                                line_new(l_counter, 1) = clone_numbering(clone_1+1)
                                line_new(l_counter, 2) = cell_numbering(host_2)
                            endif
                        elseif (clone_1 >= 0 .and. clone_2 < 0) then
                            line_new(l_counter, 2) = cell_numbering(host_2)
                            if (clone_numbering(clone_1 + 1) > 0) then
                                line_new(l_counter, 1) = clone_numbering(clone_1+1)
                            else
                                line_new(l_counter, 1) = cell_numbering(host_1)
                            endif
                        elseif (clone_1 < 0 .and. clone_2 >= 0 ) then
                            line_new(l_counter, 1) = cell_numbering(host_1)
                            if (clone_numbering(clone_2 + 1) > 0 ) then
                                line_new(l_counter, 2) = clone_numbering(clone_2+1)
                            else
                                line_new(l_counter, 2) = cell_numbering(host_2)
                            endif
                        endif
                    else
                        tot_number_lines = tot_number_lines + 1
                    endif
                endif
            endif
            if (pixel_no < num_pix) then
                new_clone_1 = clone_mask(pixel_i, new_pixel)
                new_clone_2 = clone_mask(pixel_i + 1, new_pixel)
                
                if (present(line_new)) then
                    call find_active_clone_lines(direction, pixel_i, pixel_j, min_pix, new_clone_1, new_clone_2, num_pix, area_mask, clone_mask, tot_number_lines, pixel_no, &
                                                line_new, host_1, host_2, l_counter, clone_numbering, cell_numbering)
                else
                    call find_active_clone_lines(direction, pixel_i, pixel_j, min_pix, new_clone_1, new_clone_2, num_pix, area_mask, clone_mask, tot_number_lines, pixel_no)
                endif
            endif

        elseif (direction == V_DIR) then  !! v-direction
            pixel = pixel_i + pixel_no
            do counter = pixel, pixel + num_pix - 1  !! find the common border
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
                        if (clone_1 >= 0 .and. clone_2 >= 0) then
                            if (clone_numbering(clone_1 + 1) > 0 .and. clone_numbering(clone_2 + 1) > 0) then
                                line_new(l_counter, 1) = clone_numbering(clone_1+1)
                                line_new(l_counter, 2) = clone_numbering(clone_2+1)
                            elseif (clone_numbering(clone_1 + 1) == 0 .and. clone_numbering(clone_2 + 1) > 0) then
                                line_new(l_counter, 1) = cell_numbering(host_1)
                                line_new(l_counter, 2) = clone_numbering(clone_2+1)
                            else
                                line_new(l_counter, 1) = clone_numbering(clone_1+1)
                                line_new(l_counter, 2) = cell_numbering(host_2)
                            endif
                        elseif (clone_1 >= 0 .and. clone_2 < 0) then
                            line_new(l_counter, 2) = cell_numbering(host_2)
                            if (clone_numbering(clone_1 + 1) > 0) then
                                line_new(l_counter, 1) = clone_numbering(clone_1+1)
                            else
                                line_new(l_counter, 1) = cell_numbering(host_1)
                            endif
                        elseif (clone_1 < 0 .and. clone_2 >= 0) then
                            line_new(l_counter, 1) = cell_numbering(host_1)
                            if (clone_numbering(clone_2 + 1) > 0 ) then
                                line_new(l_counter, 2) = clone_numbering(clone_2+1)
                            else
                                line_new(l_counter, 2) = cell_numbering(host_2)
                            endif
                        endif
                    else
                        tot_number_lines = tot_number_lines + 1
                    endif
                endif
            endif
            if (pixel_no < num_pix) then
                new_clone_1 = clone_mask(new_pixel, pixel_j)
                new_clone_2 = clone_mask(new_pixel, pixel_j + 1)
                if (present(line_new)) then
                    call find_active_clone_lines(direction, pixel_i, pixel_j, min_pix, new_clone_1, new_clone_2, num_pix, area_mask, clone_mask, tot_number_lines, pixel_no, &
                                                line_new, host_1, host_2, l_counter, clone_numbering, cell_numbering)
                else
                    call find_active_clone_lines(direction, pixel_i, pixel_j, min_pix, new_clone_1, new_clone_2, num_pix, area_mask, clone_mask, tot_number_lines, pixel_no)
                endif
            endif

        endif

    end subroutine find_active_clone_lines

    subroutine find_interclone_lines(clones_in_cell, l_counter, line_new, clone_numbering, clone_array)
        
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
                ! elseif (clones_in_cell(n_cell) == 2) then
            endif
        enddo

    end subroutine find_interclone_lines

end module