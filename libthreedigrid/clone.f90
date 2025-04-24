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

        !! temporary !!
        do i = 1, n_cells
            if (clone_array(i, 1) >= 0) then
                write(*,*) 'id of clone cells are ', i, clone_array(i,:)
            end if
        end do
        !! temporary !!

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
                    clone_numbering(clone_array(cell_counter, clone_counter) + 1) = 0
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
        integer*1, intent(in) :: clone_mask(:,:)             !! Identifier of the clone cells according to the area mask
        integer, intent(in) :: clones_in_cell(:)              !! Updated identifier of the clones within each host cell (does not match the clone_mask)
        integer, intent(in) :: cell_numbering(:)             !! New numbering list for quadtree cells
        integer, intent(in) :: clone_numbering(:)            !! New numbering list for clone cells
        integer, intent(inout) :: line(:,:)
        integer, intent(in) :: nodk(:)
        integer, intent(in) :: nodm(:)
        integer, intent(in) :: nodn(:)

        integer :: host_1, host_2, k_host_1, k_host_2, m_host_1, n_host_1, clone_1, clone_2
        integer :: pixel_i, pixel_j, i0, i1, j0, j1, l_counter, nl
        integer :: c, num_pix    !! temp
        integer :: n_line_new    !! temp
        integer :: pixel_no
        integer, allocatable :: cell_order(:)      ! administration of column/row no of cell wrt the direction

        write(*,*) 'size(clone_mask)', size(clone_mask,1)
        write(*,*) 'clone_mask(1441,320)', clone_mask(1441,320)
        write(*,*) 'clone_mask(1440,320)', clone_mask(1440,320)
        write(*,*) 'clone_mask(1441,321)', clone_mask(1441,321)
        write(*,*) 'clone_mask(1440,321)', clone_mask(1440,321)
        write(*,*) 'clone_mask(1440,319)', clone_mask(1440,319)
        write(*,*) 'clone_mask(1441,319)', clone_mask(1441,319)
        write(*,*) 'clone_mask(1442,319)', clone_mask(1442,319)
        write(*,*) 'clone_mask(1442,320)', clone_mask(1442,320)
        write(*,*) 'clone_mask(1442,321)', clone_mask(1442,321)

        
        n_line_new = n_line !! intent(inout)
        if (direction == U_DIR) then
            cell_order = nodm
        elseif (direction == V_DIR) then
            cell_order = nodn
        endif

        c = 0   !! temp
        do nl = start_l, end_l   !!size(line_new) n_line_u, n_line_v
            write(*,*) 'nl', nl
            if (cell_order(line(nl,1) + 1) < cell_order(line(nl,2) + 1)) then
                ! we add one (+1) to comply with fortran indexing
                host_1 = line(nl,1) + 1
                host_2 = line(nl,2) + 1
            else
                host_1 = line(nl,2) + 1
                host_2 = line(nl,1) + 1
            endif
            
            if (clones_in_cell(host_1) > 0 .or. clones_in_cell(host_2) > 0) then !! So it means that the quadtree flowline must be deleted
                c = c+1
                k_host_1 = nodk(host_1)
                k_host_2 = nodk(host_2)
                m_host_1 = nodm(host_1)
                n_host_1 = nodn(host_1)

                num_pix = lgrmin * 2**(k_host_1-1) !! temporary
                ! if (c == 1) then
                !     write(*,*) '1 - host_1', host_1, 'host_2', host_2
                !     write(*,*) 'nl', nl
                !     write(*,*) 'clones_in_cell(host_1, 1)', clones_in_cell(host_1), 'clones_in_cell(host_2, 1)',clones_in_cell(host_2)
                !     write(*,*) 'k(host_1)', k_host_1, 'k(host_2)',k_host_2
                ! endif
                if (k_host_1 <= k_host_2) then 
                    n_line_new = n_line_new - 1   !! if no refinement or refinement at the left-side, there is one flowline
                    ! if (c==1) then
                        write(*,*) '1 line is eliminated, n_line_new is ', n_line_new
                    ! endif
                else
                    n_line_new = n_line_new - 2   !! with refinement at the right-side, there are two flowlines
                    ! if (c==1) then
                        write(*,*) '2 lines are eliminated'
                    ! endif
                endif
                call get_pix_corners(k_host_1, m_host_1, n_host_1, lgrmin, i0, i1, j0, j1)
                if (direction == U_DIR) then
                    pixel_i = i1
                    pixel_j = size(clone_mask,1)-j0-(num_pix-2)
                    clone_1 = clone_mask(pixel_j,i1)
                    clone_2 = clone_mask(pixel_j,i1 + 1)
                    ! if (c == 1) then
                        write(*,*) 'clone_1', clone_1, 'clone_2', clone_2
                        write(*,*) 'i0=', i0, 'i1', i1, 'j0=', j0, 'j1', j1 
                        write(*,*) 'pixel_j', pixel_j
                    ! endif    
                elseif (direction == V_DIR) then
                    pixel_i = i0
                    pixel_j = size(clone_mask,1)-j1+num_pix
                    clone_1 = clone_mask(size(clone_mask,1)-j1, i0)
                    clone_2 = clone_mask(size(clone_mask,1)-j1 + 1,i0)    
                endif
                pixel_no = 0
                call find_active_clone_lines(direction, pixel_i, pixel_j, lgrmin, clone_1, clone_2, k_host_1, area_mask, clone_mask, n_line_new, pixel_no)   !! Only counting the lines         
                ! if (c == 1) then
                !     call find_active_clone_lines(direction, pixel_i, pixel_j, lgrmin, clone_1, clone_2, k_host_1, area_mask, clone_mask, n_line_new)   !! Only counting the lines         
                ! endif
            end if
        end do
        n_line = n_line_new

    end subroutine count_clone_lines

    subroutine find_clone_lines(direction, start_l, end_l, n_line_new, lgrmin, area_mask, clone_mask, clones_in_cell, cell_numbering, clone_numbering, line, nodk, nodm, nodn, l_counter, line_new)

        use parameters, only: U_DIR, V_DIR
        integer, intent(in) :: direction
        integer, intent(in) :: start_l
        integer, intent(in) :: end_l
        integer, intent(inout) :: n_line_new
        integer, intent(in) :: lgrmin
        integer*1, intent(in) :: area_mask(:,:)
        integer*1, intent(in) :: clone_mask(:,:)             !! Identifier of the clone cells according to the area mask
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
        integer :: pixel_i, pixel_j, i0, i1, j0, j1, nl, pixel_no
        integer, allocatable :: cell_order(:)      ! administration of column/row no of cell wrt the direction


        if (direction == U_DIR) then
            cell_order = nodm
        elseif (direction == V_DIR) then
            cell_order = nodn
        endif
    
        do nl = start_l, end_l
            if (cell_order(line(nl,1)) < cell_order(line(nl,2))) then
                host_1 = line(nl,1) + 1
                host_2 = line(nl,2) + 1
            else
                host_1 = line(nl,2) + 1
                host_2 = line(nl,1) + 1
            endif
            if (clones_in_cell(host_1) > 0 .or. clones_in_cell(host_2) > 0) then
                k_host_1 = nodk(host_1)
                k_host_2 = nodk(host_2)
                m_host_1 = nodm(host_1)
                n_host_1 = nodn(host_1)
                call get_pix_corners(k_host_1, m_host_1, n_host_1, lgrmin, i0, i1, j0, j1)
                if (direction == U_DIR) then
                    pixel_i = i1
                    pixel_j = j0
                    clone_1 = clone_mask(i1,j0)
                    clone_2 = clone_mask(i1 + 1, j0)
                elseif (direction == V_DIR) then
                    pixel_i = i0
                    pixel_j = j1
                    clone_1 = clone_mask(i0, j1)
                    clone_2 = clone_mask(i0, j1 + 1)    
                endif
                pixel_no = 0
                call find_active_clone_lines(direction, pixel_i, pixel_j, lgrmin, clone_1, clone_2, k_host_1, area_mask, clone_mask, n_line_new, pixel_no, &
                                            line_new, host_1, host_2, l_counter, clone_numbering, cell_numbering)   !! Now rewiring the line administration in u-dir
            else
                l_counter = l_counter + 1
                line_new(l_counter, 1) = cell_numbering(line(nl,1))
                line_new(l_counter, 2) = cell_numbering(line(nl,2))
            end if
        end do
    end subroutine find_clone_lines
    
    recursive subroutine find_active_clone_lines(direction, pixel_i, pixel_j, min_pix, clone_1, clone_2, k, area_mask, clone_mask, tot_number_lines, pixel_no, &
                                                line_new, host_1, host_2, l_counter, clone_numbering, cell_numbering)
        
        use parameters, only: U_DIR, V_DIR
        integer, intent(in) :: direction
        integer, intent(in) :: pixel_i
        integer, intent(in) :: pixel_j
        integer, intent(in) :: min_pix
        integer, intent(in) :: clone_1
        integer, intent(in) :: clone_2
        integer, intent(in) :: k
        integer*1, intent(in) :: area_mask(:,:)
        integer*1, intent(in) :: clone_mask(:,:)
        integer, intent(inout) :: tot_number_lines
        integer, intent(inout) :: pixel_no
        integer, intent(inout), optional :: line_new(:,:)
        integer, intent(in), optional :: host_1
        integer, intent(in), optional :: host_2
        integer, intent(inout), optional :: l_counter
        integer, intent(in), optional :: clone_numbering(:)
        integer, intent(in), optional :: cell_numbering(:)

        integer :: counter, num_pix, pixel, new_pixel
        integer :: new_clone_1, new_clone_2
        integer :: new_pixel1, new_pixel2 !!temp

        integer :: temp

        temp = 0
        num_pix = min_pix * 2**(k-1)
        ! write(*,*) 'min_pix', min_pix
        ! write(*,*) 'num_pix', num_pix
        ! pixel_no = 0
        if (direction == U_DIR) then  !! u-direction
            pixel = pixel_j + pixel_no
            do counter = pixel, pixel + num_pix - 1  !! find the common border
                write(*,*) 'counter', counter, 'pixel_i', pixel_i, 'pixel_j', pixel_j
                ! if (clone_mask(pixel_i, counter) == clone_1 .and. clone_mask(pixel_i+1, counter) == clone_2) then
                if (clone_mask( counter, pixel_i) == clone_1 .and. clone_mask( counter, pixel_i+1) == clone_2) then
                    pixel_no = pixel_no + 1
                    write(*,*) 'pixel_no', pixel_no
                    ! temp = temp + 1
                    ! if (temp==1)then
                    !     write(*,*) 'cloneline is here'
                    ! endif
                else
                    exit
                endif
            enddo
            new_pixel1 = pixel_j + pixel_no !! I can also use counter instead of new_pixel
            new_pixel2 = size(clone_mask,1) - new_pixel1 + 2
            write(*,*) 'pixel no area-mask version:', new_pixel2
            write(*,*) 'pixel no clone-mask version:', new_pixel1
            write(*,*) 'pixel no start for area mask:', size(clone_mask,1)-pixel+1
            ! if (any(minval(area_mask(pixel_i:pixel_i+1,pixel_j:new_pixel2-1), 1) > 0)) then  !! if there is a minimum of a pair of pixel with data, make the flowline
            if (any(minval(area_mask(pixel_i:pixel_i+1,new_pixel2:(size(clone_mask,1)-pixel+1)), 1) > 0)) then  !! if there is a minimum of a pair of pixel with data, make the flowline
                write(*,*) 'yes, create a new line'
                if (present(line_new)) then
                    if (clone_1 > 0 .and. clone_2 > 0) then
                        l_counter = l_counter + 1
                        line_new(l_counter, 1) = clone_numbering(clone_1+1)
                        line_new(l_counter, 2) = clone_numbering(clone_2+1)
                    elseif (clone_1 > 0 .and. clone_2 <= 0) then
                        l_counter = l_counter + 1
                        line_new(l_counter, 1) = clone_numbering(clone_1+1)
                        line_new(l_counter, 2) = cell_numbering(host_2)
                    elseif (clone_2 > 0 .and. clone_1 <= 0) then
                        l_counter = l_counter + 1
                        line_new(l_counter, 1) = clone_numbering(host_1)
                        line_new(l_counter, 2) = cell_numbering(clone_2+1)
                    endif
                    ! write(*,*) 'line_new_1', line_new(counter, 1)
                    ! write(*,*) 'line_new_2', line_new(counter, 2)
                else
                    tot_number_lines = tot_number_lines + 1

                    write(*,*) '1 line is created, tot_number_lines is ', tot_number_lines
                endif
            endif
            write(*,*) 'pixel_no', pixel_no, 'num_pix', num_pix
            if (pixel_no < num_pix) then
                write(*,*) 'more than 1 line is created'
                ! new_clone_1 = clone_mask(pixel_i, new_pixel1)
                ! new_clone_2 = clone_mask(pixel_i + 1, new_pixel1)
                new_clone_1 = clone_mask(new_pixel1, pixel_i)
                new_clone_2 = clone_mask( new_pixel1, pixel_i + 1)
                if (present(line_new)) then
                    call find_active_clone_lines(direction, pixel_i, pixel_j, min_pix, new_clone_1, new_clone_2, k, area_mask, clone_mask, tot_number_lines, pixel_no, &
                                                line_new, host_1, host_2, l_counter, clone_numbering, cell_numbering)
                else
                    call find_active_clone_lines(direction, pixel_i, pixel_j, min_pix, new_clone_1, new_clone_2, k, area_mask, clone_mask, tot_number_lines, pixel_no)
                endif
            endif

        elseif (direction == V_DIR) then  !! v-direction
            pixel = pixel_i
            do counter = pixel, pixel + num_pix - 1  !! find the common border
                ! if (clone_mask(counter, pixel_j) == clone_1 .and. clone_mask(counter, pixel_j+1) == clone_2) then
                if (clone_mask(pixel_j, counter) == clone_1 .and. clone_mask(pixel_j+1,counter) == clone_2) then
                    pixel_no = pixel_no + 1
                else
                    exit
                endif
            enddo
            new_pixel = pixel + pixel_no
            if (any(minval(area_mask(pixel_i:new_pixel-1,pixel_j:pixel_j+1), 2) > 0)) then
                if (present(line_new)) then
                    if (clone_1 > 0 .and. clone_2 > 0) then
                        l_counter = l_counter + 1
                        line_new(l_counter, 1) = clone_numbering(clone_1+1)
                        line_new(l_counter, 2) = clone_numbering(clone_2+1)
                    elseif (clone_1 > 0 .and. clone_2 <= 0) then
                        l_counter = l_counter + 1
                        line_new(l_counter, 1) = clone_numbering(clone_1+1)
                        line_new(l_counter, 2) = cell_numbering(host_2)
                    elseif (clone_2 > 0 .and. clone_1 <= 0) then
                        l_counter = l_counter + 1
                        line_new(l_counter, 1) = cell_numbering(host_1)
                        line_new(l_counter, 2) = clone_numbering(clone_2+1)
                    endif
                else
                    tot_number_lines = tot_number_lines + 1
                endif
            endif
            if (pixel_no < num_pix) then
                ! new_clone_1 = clone_mask(new_pixel, pixel_j)
                ! new_clone_2 = clone_mask(new_pixel, pixel_j + 1)
                new_clone_1 = clone_mask(pixel_j,new_pixel)
                new_clone_2 = clone_mask(pixel_j + 1,new_pixel)
                if (present(line_new)) then
                    call find_active_clone_lines(direction, new_pixel, pixel_j, min_pix, new_clone_1, new_clone_2, k, area_mask, clone_mask, tot_number_lines, pixel_no, &
                                                line_new, host_1, host_2, l_counter, clone_numbering, cell_numbering)
                else
                    call find_active_clone_lines(direction, new_pixel, pixel_j, min_pix, new_clone_1, new_clone_2, k, area_mask, clone_mask, tot_number_lines, pixel_no)
                endif
            endif

        endif

    end subroutine find_active_clone_lines

end module