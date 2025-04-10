module m_clone

    use m_grid_utils, only : get_pix_corners
    
    implicit none

    contains

    subroutine find_active_clone_cells(n_cells, clone_array, cell_numbering, clone_numbering)  !!, n_clone_cells, clones_in_cell
        ! ! Finding clone cells based on the area mask specifically defined for clone cells
        
        use m_grid_utils, only : get_pix_corners
        use parameters, only : CLONE_NUMBERS
        integer, intent(inout) :: n_cells                  !! total number of active cells. in: without clone cells, out: with clone cells
        integer, intent(inout) :: clone_array(:,:)         !! Identifier of the clones within each host cell. in: matches the clone_mask, out: matches the new numbering
        integer, intent(inout) :: cell_numbering(:)            !! New numbering for quadtree cells
        integer, intent(inout) :: clone_numbering(:)            !! New numbering for clone cells
        integer :: clones_in_cell(n_cells)                       !! Total number of clones in each host cell
        integer :: cell_counter, clone_counter, i, counter
        integer :: n_cells_new

        ! write(*,*) 'fortran', clone_array(107,0)
        write(*,*) 'clone_array(n_cells,:)', clone_array(n_cells,:)
        ! write(*,*) 'a(1368,479) = ', clone_mask(1368,479)
        ! write(*,*) 'max(clone_array) = ', maxval(clone_array)
        ! write(*,*) 'a(481,1368) = ', clone_mask(481,1368)
        do i = 1, n_cells
            if (clone_array(i, 1) >= 0) then
                write(*,*) 'id of clone cells are ', clone_array(i,:)
            end if
        end do
        counter = 0
        ! write(*,*) 'fortran - clone_array(0,0)', clone_array(0, 0)
        ! write(*,*) 'fortran - clone_array(0,1)', clone_array(0, 1)
        write(*,*) 'fortran - clone_array(1,1)', clone_array(1, 1)
        write(*,*) 'max(clone_array) = ', maxval(clone_array)

        ! n_clone_cells = maxval(clone_mask)
        n_cells_new = n_cells
        do cell_counter = 1, n_cells
            if (clone_array(cell_counter, 1) >= 0) then
                do clone_counter = 1, CLONE_NUMBERS
                    if (clone_array(cell_counter, clone_counter) >= 0) then
                        counter = counter + 1
                        clone_numbering(clone_array(cell_counter, clone_counter)) = counter
                        clones_in_cell(cell_counter) = clones_in_cell(cell_counter) + 1
                    end if
                end do
            else
                counter = counter + 1
                cell_numbering(cell_counter) = counter
            end if
            if (clones_in_cell(cell_counter) > 0) then
                n_cells_new = n_cells_new + clones_in_cell(cell_counter) - 1
            end if
        end do

        n_cells = n_cells_new

    end subroutine find_active_clone_cells

    subroutine find_clone_lines(n_line_u, n_line_v, lgrmin, area_mask, clone_mask, clone_array, cell_numbering, clone_numbering, line, nodk, nodm, nodn)

        integer, intent(inout) :: n_line_u
        integer, intent(inout) :: n_line_v
        integer, intent(inout) :: lgrmin
        integer*1, intent(in) :: area_mask(:,:)
        integer*1, intent(in) :: clone_mask(:,:)             !! Identifier of the clone cells according to the area mask
        integer, intent(in) :: clone_array(:,:)              !! Updated identifier of the clones within each host cell (does not match the clone_mask)
        integer, intent(in) :: cell_numbering(:)             !! New numbering list for quadtree cells
        integer, intent(in) :: clone_numbering(:)            !! New numbering list for clone cells
        integer, intent(inout) :: line(:,:)
        integer, intent(in) :: nodk(:)
        integer, intent(in) :: nodm(:)
        integer, intent(in) :: nodn(:)

        integer :: new_line_u, new_line_v, n_cloneline_u, n_cloneline_v, nl
        integer :: host_1, host_2, k_host_1, k_host_2, m_host, n_host, clone_1, clone_2
        integer :: i0, i1, j0, j1, l_counter
        integer, allocatable :: line_new(:,:)      ! administration of clone lines

        new_line_u = n_line_u
        n_cloneline_u = 0
        do nl = 1, n_line_u   !!size(line_new) n_line_u, n_line_v
            if (nodm(line(nl,1)) < nodm(line(nl,2))) then
                host_1 = line(nl,1)
                host_2 = line(nl,2)
            else
                host_1 = line(nl,2)
                host_2 = line(nl,1)
            endif
            k_host_1 = nodk(host_1)
            k_host_2 = nodk(host_2)
            m_host = nodm(host_1)
            n_host = nodn(host_1)
            if (clone_array(host_1, 1) > 0 .or. clone_array(host_2, 1) > 0) then !! So it means that the quadtree flowline must be deleted
                write(*,*) 'nl', nl
                write(*,*) 'clone_array(host_1, 1)', clone_array(host_1, 1), 'clone_array(host_2, 1)',clone_array(host_2, 1)
                if (k_host_1 <= k_host_2) then 
                    new_line_u = new_line_u - 1   !! if no refinement, there is only one flowline
                else
                    new_line_u = new_line_u - 2   !! with refinement, there are two flowlines
                endif
                call get_pix_corners(k_host_1, m_host, n_host, lgrmin, i0, i1, j0, j1)
                clone_1 = clone_mask(i1,j0)
                clone_2 = clone_mask(i1 + 1, j0)
                
                call find_active_clone_lines(1, i1, j0, lgrmin, clone_1, clone_2, k_host_1, area_mask, clone_mask, new_line_u)   !! Only counting the lines         
            end if
        end do

        new_line_v = n_line_v
        n_cloneline_v = 0
        do nl = n_line_u + 1, n_line_u + n_line_v
            if (nodn(line(nl,1)) < nodn(line(nl,2))) then
                host_1 = line(nl,1)
                host_2 = line(nl,2)
            else
                host_1 = line(nl,2)
                host_2 = line(nl,1)
            endif
            k_host_1 = nodk(host_1)
            k_host_2 = nodk(host_2)
            m_host = nodm(host_1)
            n_host = nodn(host_1)
            if (clone_array(host_1, 1) > 0 .or. clone_array(host_2, 1) > 0) then
                if (k_host_1 == k_host_2) then !! So it means that the quadtree flowline must be deleted
                    new_line_v = new_line_v - 1   !! if no refinement, there is only one flowline
                else
                    new_line_v = new_line_v - 2   !! with refinement, there are two flowlines
                endif
                call get_pix_corners(k_host_1, m_host, n_host, lgrmin, i0, i1, j0, j1)
                clone_1 = clone_mask(i0, j1)
                clone_2 = clone_mask(i0, j1 + 1)

                call find_active_clone_lines(2, i1, j0, lgrmin, clone_1, clone_2, k_host_1, area_mask, clone_mask, new_line_v)    !! Only counting the lines         
            end if
        end do

        allocate(line_new(new_line_u + new_line_v, 2)) !! this includes all lines with new numbering for the cells
        l_counter = 0
        do nl = 1, n_line_u
            if (nodm(line(nl,1)) < nodm(line(nl,2))) then
                host_1 = line(nl,1)
                host_2 = line(nl,2)
            else
                host_1 = line(nl,2)
                host_2 = line(nl,1)
            endif
            k_host_1 = nodk(host_1)
            k_host_2 = nodk(host_2)
            m_host = nodm(host_1)
            n_host = nodn(host_1)
            if (clone_array(host_1, 1) > 0 .or. clone_array(host_2, 1) > 0) then
                call get_pix_corners(k_host_1, m_host, n_host, lgrmin, i0, i1, j0, j1)
                clone_1 = clone_mask(i1,j0)
                clone_2 = clone_mask(i1 + 1, j0)

                call find_active_clone_lines(1, i1, j0, lgrmin, clone_1, clone_2, k_host_1, area_mask, clone_mask, new_line_u, &
                                            line_new, host_1, host_2, l_counter, clone_numbering, cell_numbering)   !! Now rewiring the line administration in u-dir
            else
                l_counter = l_counter + 1
                line_new(l_counter, 1) = cell_numbering(line(nl,1))
                line_new(l_counter, 2) = cell_numbering(line(nl,2))
            end if
        end do

        do nl = n_line_u + 1, n_line_u + n_line_v
            if (nodn(line(nl,1)) < nodn(line(nl,2))) then
                host_1 = line(nl,1)
                host_2 = line(nl,2)
            else
                host_1 = line(nl,2)
                host_2 = line(nl,1)
            endif
            k_host_1 = nodk(host_1)
            k_host_2 = nodk(host_2)
            m_host = nodm(host_1)
            n_host = nodn(host_1)
            if (clone_array(host_1, 1) > 0 .or. clone_array(host_2, 1) > 0) then
                call get_pix_corners(k_host_1, m_host, n_host, lgrmin, i0, i1, j0, j1)
                clone_1 = clone_mask(i0, j1)
                clone_2 = clone_mask(i0, j1 + 1)

                call find_active_clone_lines(2, i1, j0, lgrmin, clone_1, clone_2, k_host_1, area_mask, clone_mask, new_line_v, &
                                            line_new, host_1, host_2, l_counter, clone_numbering, cell_numbering)    !! Now rewiring the line administration in v-dir
            else
                l_counter = l_counter + 1
                line_new(l_counter, 1) = cell_numbering(line(nl,1))
                line_new(l_counter, 2) = cell_numbering(line(nl,2))
            end if
        end do

        write(*,*) 'line_new', line_new(:,:)
        n_line_u = new_line_u
        n_line_v = new_line_v
        line = line_new

    end subroutine find_clone_lines

    recursive subroutine find_active_clone_lines(direction, ii, jj, min_pix, clone_1, clone_2, k, area_mask, clone_mask, tot_number_lines, &
                                                line_new, host_1, host_2, l_new, clone_numbering, cell_numbering)
                    
        integer, intent(in) :: direction
        integer, intent(in) :: ii
        integer, intent(in) :: jj
        integer, intent(in) :: min_pix
        integer, intent(in) :: clone_1
        integer, intent(in) :: clone_2
        integer, intent(in) :: k
        integer*1, intent(in) :: area_mask(:,:)
        integer*1, intent(in) :: clone_mask(:,:)
        integer, intent(inout) :: tot_number_lines
        integer, intent(inout), optional :: line_new(:,:)
        integer, intent(in), optional :: host_1
        integer, intent(in), optional :: host_2
        integer, intent(inout), optional :: l_new
        integer, intent(in), optional :: clone_numbering(:)
        integer, intent(in), optional :: cell_numbering(:)

        integer :: pixel_no, counter, num_pix, pixel, new_pixel
        integer :: new_clone_1, new_clone_2

        num_pix = min_pix * 2**(k-1)
        write(*,*) 'min_pix', min_pix
        write(*,*) 'num_pix', num_pix
        pixel_no = 0
        if (direction == 1) then  !! u-direction
            pixel = jj
            do counter = pixel, pixel + num_pix - 1
                if (clone_mask(ii,counter) == clone_1 .and. clone_mask(ii+1,counter) == clone_2) then
                    pixel_no = pixel_no + 1
                else
                    exit
                endif
            enddo
            new_pixel = pixel + pixel_no !! I can also use counter instead of new_pixel
            if (any(minval(area_mask(ii:ii+1,jj:new_pixel-1), 1) > 0)) then
                if (present(line_new)) then
                    if (clone_1 > 0 .and. clone_2 > 0) then
                        l_new = l_new + 1
                        line_new(l_new, 1) = clone_numbering(clone_1)
                        line_new(l_new, 2) = clone_numbering(clone_2)
                    elseif (clone_1 > 0 .and. clone_2 <= 0) then
                        l_new = l_new + 1
                        line_new(l_new, 1) = clone_numbering(clone_1)
                        line_new(l_new, 2) = cell_numbering(host_2)
                    elseif (clone_2 > 0 .and. clone_1 <= 0) then
                        l_new = l_new + 1
                        line_new(l_new, 1) = clone_numbering(host_1)
                        line_new(l_new, 2) = cell_numbering(clone_2)
                    endif
                    write(*,*) 'line_new_1', line_new(counter, 1)
                    write(*,*) 'line_new_2', line_new(counter, 2)
                else
                    tot_number_lines = tot_number_lines + 1
                    write(*,*) 'tot_number_lines', tot_number_lines
                endif
            endif
            if (pixel_no < num_pix) then
                new_clone_1 = clone_mask(ii, new_pixel)
                new_clone_2 = clone_mask(ii + 1, new_pixel)
                if (present(line_new)) then
                    call find_active_clone_lines(direction, ii, new_pixel, min_pix, new_clone_1, new_clone_2, k, area_mask, clone_mask, tot_number_lines, &
                                                line_new, host_1, host_2, l_new, clone_numbering, cell_numbering)
                else
                    call find_active_clone_lines(direction, ii, new_pixel, min_pix, new_clone_1, new_clone_2, k, area_mask, clone_mask, tot_number_lines)
                endif
            endif

        elseif (direction == 2) then  !! v-direction
            pixel = ii
            do counter = pixel, pixel + num_pix - 1
                if (clone_mask(counter,jj) == clone_1 .and. clone_mask(counter,jj+1) == clone_2) then
                    pixel_no = pixel_no + 1
                else
                    exit
                endif
            enddo
            new_pixel = pixel + pixel_no !! I can also use counter instead of new_pixel
            if (any(minval(area_mask(ii:new_pixel-1,jj:jj+1), 2) > 0)) then
                if (present(line_new)) then
                    if (clone_1 > 0 .and. clone_2 > 0) then
                        l_new = l_new + 1
                        line_new(l_new, 1) = clone_numbering(clone_1)
                        line_new(l_new, 2) = clone_numbering(clone_2)
                    elseif (clone_1 > 0 .and. clone_2 <= 0) then
                        l_new = l_new + 1
                        line_new(l_new, 1) = clone_numbering(clone_1)
                        line_new(l_new, 2) = cell_numbering(host_2)
                    elseif (clone_2 > 0 .and. clone_1 <= 0) then
                        l_new = l_new + 1
                        line_new(l_new, 1) = cell_numbering(host_1)
                        line_new(l_new, 2) = clone_numbering(clone_2)
                    endif
                else
                    tot_number_lines = tot_number_lines + 1
                endif
            endif
            if (pixel_no < num_pix) then
                new_clone_1 = clone_mask(new_pixel, jj)
                new_clone_2 = clone_mask(new_pixel, jj + 1)
                if (present(line_new)) then
                    call find_active_clone_lines(direction, new_pixel, jj, min_pix, new_clone_1, new_clone_2, k, area_mask, clone_mask, tot_number_lines, &
                                                line_new, host_1, host_2, l_new, clone_numbering, cell_numbering)
                else
                    call find_active_clone_lines(direction, new_pixel, jj, min_pix, new_clone_1, new_clone_2, k, area_mask, clone_mask, tot_number_lines)
                endif
            endif

        endif

    end subroutine find_active_clone_lines

end module