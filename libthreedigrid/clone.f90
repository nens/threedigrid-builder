module m_clone

    implicit none

    contains

    subroutine find_active_clone_cells(cloned, area_mask_clone, n_cells, n_clone_cells)  !!, n_clone_cells, clones_in_cell
        ! ! Finding clone cells based on the area mask specifically defined for clone cells
        
        use m_grid_utils, only : get_pix_corners
        use parameters, only : CLONE_NUMBERS
        integer, intent(inout) :: n_cells                  !! total number of active cells
        integer, intent(in) :: cloned(:,:)              !! Identifier of the clones within each host cell
        integer*1, intent(in) :: area_mask_clone(:,:)   !! Identifier of the clone cells according to the area mask
        integer, intent(inout) :: n_clone_cells     !! Total number of clone cells
        ! integer, intent(inout) :: clones_in_cell(:)         !! number of clones in each host cell
        integer :: cell_counter, clone_counter, i
        integer :: n_cells_new
        integer :: clones_in_cell(n_cells)

        ! write(*,*) 'fortran', cloned(107,0)
        ! write(*,*) cloned(0,:)
        ! write(*,*) 'a(1368,479) = ', area_mask_clone(1368,479)
        ! write(*,*) 'max(cloned) = ', maxval(cloned)
        ! write(*,*) 'a(481,1368) = ', area_mask_clone(481,1368)
        ! do i = 1, n_cells
        !     if (cloned(1, i) >= 0) then
        !         write(*,*) 'id of clone cells are ', cloned(:,i)
        !     end if
        ! end do
        n_clone_cells = 0
        n_cells_new = n_cells
        do cell_counter = 0, n_cells-1
            if (cloned(cell_counter, 1) >= 0) then
                do clone_counter = 1, CLONE_NUMBERS
                    if (cloned(cell_counter, clone_counter) >= 0) then
                        n_clone_cells = n_clone_cells + 1
                        clones_in_cell(cell_counter) = clones_in_cell(cell_counter) + 1
                    end if
                end do
            end if
            if (clones_in_cell(cell_counter) > 0) then
                n_cells_new = n_cells_new + clones_in_cell(cell_counter) - 1
            end if
        end do

        n_cells = n_cells_new

    end subroutine find_active_clone_cells

    subroutine find_clone_lines(n_line_u, n_line_v,lgrmin,area_mask,area_mask_clone,cloned,line,nodk,nodm,nodn)

        integer, intent(inout) :: n_line_u
        integer, intent(inout) :: n_line_v
        integer, intent(inout) :: lgrmin
        integer*1, intent(in) :: area_mask(:,:)
        integer*1, intent(in) :: area_mask_clone(:,:)   !! Identifier of the clone cells according to the area mask
        integer, intent(in) :: cloned(:,:)              !! Identifier of the clones within each host cell
        integer, intent(in) :: line(:,:)
        integer, intent(in) :: nodk(:)
        integer, intent(in) :: nodm(:)
        integer, intent(in) :: nodn(:)

        integer :: new_line_u, new_line_v, n_cloneline_u, n_cloneline_v, nl
        integer :: host_1, host_2, k_host_1, k_host_2, m_host, n_host, clone_1, clone_2
        integer :: i0, i1, j0, j1
        integer, allocatable :: clone_line(:,:)      ! Wet perimiter area

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
            if (cloned(host_1, 1) > 0 .or. cloned(host_2, 1) > 0) then
                write(*,*) 'nl', nl
                write(*,*) 'cloned(host_1, 1)', cloned(host_1, 1), 'cloned(host_2, 1)',cloned(host_2, 1)
                if (k_host_1 <= k_host_2) then !! So it means that the quadtree flowline must be deleted
                    new_line_u = new_line_u - 1   !! if no refinement, there is only one flowline
                else
                    new_line_u = new_line_u - 2   !! with refinement, there are two flowlines
                endif
                call get_pix_corners(k_host_1, m_host, n_host, lgrmin, i0, i1, j0, j1)
                clone_1 = area_mask_clone(i1,j0)
                clone_2 = area_mask_clone(i1 + 1, j0)
                
                call find_active_clone_lines(1, i1, j0, lgrmin, clone_1, clone_2, k_host_1, area_mask, area_mask_clone, new_line_u, n_cloneline_u)   !! Only counting the clone lines         
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
            if (cloned(host_1, 1) > 0 .or. cloned(host_2, 1) > 0) then
                if (k_host_1 == k_host_2) then !! So it means that the quadtree flowline must be deleted
                    new_line_v = new_line_v - 1   !! if no refinement, there is only one flowline
                else
                    new_line_v = new_line_v - 2   !! with refinement, there are two flowlines
                endif
                call get_pix_corners(k_host_1, m_host, n_host, lgrmin, i0, i1, j0, j1)
                clone_1 = area_mask_clone(i0, j1)
                clone_2 = area_mask_clone(i0, j1 + 1)

                call find_active_clone_lines(2, i1, j0, lgrmin, clone_1, clone_2, k_host_1, area_mask, area_mask_clone, new_line_v, n_cloneline_v)    !! Only counting the clone lines         
            end if
        end do

        ! allocate(clone_line(new_line_u + new_line_v, 2))  !! This will be a global parameter and an argument
        allocate(clone_line(n_cloneline_u + n_cloneline_v, 2))  !! This will be a global parameter and an argument
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
            if (cloned(host_1, 1) > 0 .or. cloned(host_2, 1) > 0) then
                call get_pix_corners(k_host_1, m_host, n_host, lgrmin, i0, i1, j0, j1)
                clone_1 = area_mask_clone(i1,j0)
                clone_2 = area_mask_clone(i1 + 1, j0)

                call find_active_clone_lines(1, i1, j0, lgrmin, clone_1, clone_2, k_host_1, area_mask, area_mask_clone, new_line_u, n_cloneline_u, clone_line, host_1, host_2, nl)   !! Only counting the clone lines         
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
            if (cloned(host_1, 1) > 0 .or. cloned(host_2, 1) > 0) then
                call get_pix_corners(k_host_1, m_host, n_host, lgrmin, i0, i1, j0, j1)
                clone_1 = area_mask_clone(i0, j1)
                clone_2 = area_mask_clone(i0, j1 + 1)

                call find_active_clone_lines(2, i1, j0, lgrmin, clone_1, clone_2, k_host_1, area_mask, area_mask_clone, new_line_v, n_cloneline_v, clone_line, host_1, host_2, nl)    !! Only counting the clone lines         
            end if
        end do

        write(*,*) 'clone_line', clone_line(:,:)
        n_line_u = new_line_u
        n_line_v = new_line_v

    end subroutine find_clone_lines

    recursive subroutine find_active_clone_lines(direction, ii, jj, min_pix, clone_1, clone_2, k, area_mask, area_mask_clone, tot_number_lines, number_clonelines, clone_line, host_1, host_2, l)
                    
        integer, intent(in) :: direction
        integer, intent(in) :: ii
        integer, intent(in) :: jj
        integer, intent(in) :: min_pix
        integer, intent(in) :: clone_1
        integer, intent(in) :: clone_2
        integer, intent(in) :: k
        integer*1, intent(in) :: area_mask(:,:)
        integer*1, intent(in) :: area_mask_clone(:,:)
        integer, intent(inout) :: tot_number_lines
        integer, intent(inout) :: number_clonelines
        integer, intent(inout), optional :: clone_line(:,:)
        integer, intent(in), optional :: host_1
        integer, intent(in), optional :: host_2
        integer, intent(in), optional :: l

        integer :: pixel_no, counter, num_pix, pixel, new_pixel
        integer :: new_clone_1, new_clone_2

        num_pix = min_pix * 2**(k-1)
        write(*,*) 'min_pix', min_pix
        write(*,*) 'num_pix', num_pix
        pixel_no = 0
        if (direction == 1) then  !! u-direction
            pixel = jj
            do counter = pixel, pixel + num_pix - 1
                if (area_mask_clone(ii,counter) == clone_1 .and. area_mask_clone(ii+1,counter) == clone_2) then
                    pixel_no = pixel_no + 1
                else
                    exit
                endif
            enddo
            new_pixel = pixel + pixel_no !! I can also use counter instead of new_pixel
            if (any(minval(area_mask(ii:ii+1,jj:new_pixel-1), 1) > 0)) then
                if (present(clone_line)) then
                    if (clone_1 > 0 .and. clone_2 > 0) then
                        clone_line(l, 1) = clone_1
                        clone_line(l, 2) = clone_2
                    elseif (clone_1 > 0 .and. clone_2 <= 0) then
                        clone_line(l, 1) = clone_1
                        clone_line(l, 2) = host_2
                    elseif (clone_2 > 0 .and. clone_1 <= 0) then
                        clone_line(l, 1) = host_1
                        clone_line(l, 2) = clone_2
                    endif
                    write(*,*) 'clone_line_1', clone_line(counter, 1)
                    write(*,*) 'clone_line_2', clone_line(counter, 2)
                else
                    tot_number_lines = tot_number_lines + 1
                    number_clonelines = number_clonelines + 1
                    write(*,*) 'number_clonelines', number_clonelines
                endif
            endif
            if (pixel_no < num_pix) then
                new_clone_1 = area_mask_clone(ii, new_pixel)
                new_clone_2 = area_mask_clone(ii + 1, new_pixel)
                if (present(clone_line)) then
                    call find_active_clone_lines(direction, ii, new_pixel, min_pix, new_clone_1, new_clone_2, k, area_mask, area_mask_clone, tot_number_lines, number_clonelines, clone_line, host_1, host_2)
                else
                    call find_active_clone_lines(direction, ii, new_pixel, min_pix, new_clone_1, new_clone_2, k, area_mask, area_mask_clone, tot_number_lines, number_clonelines)
                endif
            endif

        elseif (direction == 2) then  !! v-direction
            pixel = ii
            do counter = pixel, pixel + num_pix - 1
                if (area_mask_clone(counter,jj) == clone_1 .and. area_mask_clone(counter,jj+1) == clone_2) then
                    pixel_no = pixel_no + 1
                else
                    exit
                endif
            enddo
            new_pixel = pixel + pixel_no !! I can also use counter instead of new_pixel
            if (any(minval(area_mask(ii:new_pixel-1,jj:jj+1), 2) > 0)) then
                if (present(clone_line)) then
                    if (clone_1 > 0 .and. clone_2 > 0) then
                        clone_line(l, 1) = clone_1
                        clone_line(l, 2) = clone_2
                    elseif (clone_1 > 0 .and. clone_2 <= 0) then
                        clone_line(l, 1) = clone_1
                        clone_line(l, 2) = host_2
                    elseif (clone_2 > 0 .and. clone_1 <= 0) then
                        clone_line(l, 1) = host_1
                        clone_line(l, 2) = clone_2
                    endif
                else
                    tot_number_lines = tot_number_lines + 1
                    number_clonelines = number_clonelines + 1
                endif
            endif
            if (pixel_no < num_pix) then
                new_clone_1 = area_mask_clone(new_pixel, jj)
                new_clone_2 = area_mask_clone(new_pixel, jj + 1)
                if (present(clone_line)) then
                    call find_active_clone_lines(direction,new_pixel, jj, min_pix, new_clone_1, new_clone_2, k, area_mask, area_mask_clone, tot_number_lines, number_clonelines, clone_line, host_1, host_2)
                else
                    call find_active_clone_lines(direction,new_pixel, jj, min_pix, new_clone_1, new_clone_2, k, area_mask, area_mask_clone, tot_number_lines, number_clonelines)
                endif
            endif

        endif

    end subroutine find_active_clone_lines

end module