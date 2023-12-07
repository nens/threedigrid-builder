module m_quadtree

    use parameters, only : NODATA

    implicit none

    contains

    subroutine make_quadtree(kmax, mmax, nmax, lgrmin, use_2d_flow, area_mask, lg, quad_idx,&
        n_cells, n_line_u, n_line_v)
    !!! Entry point for creating quadtree by setting lg array and then finding number of active cells and lines.
        use parameters, only : NODATA

        integer, intent(in) :: kmax ! Maximum refinement levels
        integer, intent(in) :: mmax(:) ! X Dimension of each refinement level
        integer, intent(in) :: nmax(:) ! Y Dimension of each refinement level
        integer, intent(in) :: lgrmin ! Number of pixels in cell of smallest refinement level
        integer, intent(in) :: use_2d_flow  ! Whether to add flowlines
        integer*1, intent(in) :: area_mask(:, :) ! Array with active pixels of model.
        integer, intent(inout) :: lg(:, :) ! Array with all refinement levels.
        integer, intent(inout) :: quad_idx(:, :) ! Array with idx of cell at lg refinement locations
        integer, intent(inout) :: n_cells ! counter for active cells
        integer, intent(inout) :: n_line_u ! counter for active u lines
        integer, intent(inout) :: n_line_v ! counter for active v lines
        integer :: k
        integer :: m, n

        write(*,*) '** INFO: Start making quadtree.'
        do m=1, mmax(kmax)
            do n=1, nmax(kmax)
                call divide(kmax, m, n, lg)
            enddo
        enddo
        call balance_quadtree(kmax, mmax, nmax, lg)
        call find_active_2d_comp_cells(&
            kmax, mmax, nmax, lgrmin, use_2d_flow > 0, lg, area_mask, quad_idx, n_cells, n_line_u, n_line_v&
        )
        write(*,*) '** INFO: Done making quadtree.'

    end subroutine make_quadtree

    recursive subroutine divide(k, m, n, lg) !ip, jp, 
    !!! Recursive subroutine to set correct refinement levels on lg refinement array.
        use m_grid_utils, only : get_lg_corners
        use parameters, only : NODATA

        integer, intent(in) :: k
        integer, intent(in) :: m
        integer, intent(in) :: n
        integer, intent(inout) :: lg(:,:)
        integer :: k1
        integer :: mn(4)
        integer :: m_n, n_n

        mn = get_lg_corners(k, m, n)
        if(any(lg(mn(1):mn(3),mn(2):mn(4)) < k).and.k>1) then
            k1 = k-1
            m_n=2*m-1
            n_n=2*n-1
            lg(mn(1):mn(3),mn(2):mn(4)) = min(lg(mn(1):mn(3),mn(2):mn(4)), k1)
            call divide(k1, m_n, n_n, lg)
            call divide(k1, m_n+1, n_n, lg)
            call divide(k1, m_n, n_n+1, lg)
            call divide(k1, m_n+1, n_n+1, lg)
        endif

    end subroutine divide

    subroutine balance_quadtree(kmax, mmax, nmax, lg)
    !!! Balancing out refinement levels, so that neighbouring cells are never more than 1 refinement level apart.
        integer, intent(in) :: kmax
        integer, intent(in) :: mmax(:)
        integer, intent(in) :: nmax(:)
        integer, intent(inout) :: lg(:,:)
        integer :: k, k1
        integer :: m, n, m1, n1

        do k=1,kmax-1
            do n=1, nmax(1)
                do m=1,mmax(1)-1
                    if(lg(m,n) == k .and. lg(m+1,n)>k+1) then
                        lg(m+1,n)=k+1
                    endif
                    if(lg(m,n) > k+1 .and. lg(m+1,n)==k) then
                        lg(m,n)=k+1
                    endif
                enddo
            enddo
            do n=1,nmax(1)-1
                do m=1, mmax(1)
                    if (lg(m,n)==k.and.lg(m,n+1)>k+1) then
                        lg(m,n+1)=k+1
                    endif
                    if (lg(m,n)>k+1.and.lg(m,n+1)==k) then
                        lg(m,n)=k+1
                    endif
                enddo
            enddo

            
            do n=1, nmax(kmax)
                do m=1, mmax(kmax)
                    call divide(kmax, m, n, lg)
                enddo
            enddo
        enddo

    end subroutine balance_quadtree

    subroutine find_active_2d_comp_cells(&
        kmax, mmax, nmax, lgrmin, use_2d_flow, lg, area_mask, quad_idx, n_cells, n_line_u, n_line_v&
    )
    !!! Counting active cells and lines based on area_mask of active pixels.
        use m_grid_utils, only : get_lg_corners, get_pix_corners, crop_pix_coords_to_raster, pad_area_mask
        use m_cells, only : set_2d_computational_lines
        use parameters, only : NODATA

        integer, intent(in) :: kmax
        integer, intent(in) :: mmax(:)
        integer, intent(in) :: nmax(:)
        integer, intent(in) :: lgrmin
        logical, intent(in) :: use_2d_flow
        integer, intent(inout) :: lg(:,:)
        integer*1, intent(in) :: area_mask(:,:)
        integer, intent(inout) :: quad_idx(:,:)
        integer, intent(inout) :: n_line_u
        integer, intent(inout) :: n_line_v
        integer*1, allocatable:: area_mask_padded(:, :)
        integer :: k
        integer :: m,n
        integer :: mn(4)
        integer :: i0, i1, j0, j1, i2, i3, j2, j3
        integer :: n_cells
        
        
        n_cells = 0
        n_line_u = 0
        n_line_v = 0
        quad_idx = 0
        call get_pix_corners(kmax, mmax(kmax), nmax(kmax), lgrmin, i0, i1, j0, j1)
        area_mask_padded = pad_area_mask(area_mask, i0, i1, j0, j1) 
        do k=kmax,1,-1
            do m=1,mmax(k)
                do n=1,nmax(k)
                    call get_pix_corners(k, m, n, lgrmin, i0, i1, j0, j1)
                    mn = get_lg_corners(k, m, n)
                    i1 = min(i1, size(area_mask, 1))
                    j1 = min(j1, size(area_mask, 2))
                    if (all(lg(mn(1):mn(3),mn(2):mn(4)) == k)) then !! TODO: CHECK OF MODEL AREA CHECK IS NECESSARY???
                        if (all(area_mask_padded(i0:i1, j0:j1) == 0)) then
                            lg(mn(1):mn(3),mn(2):mn(4)) = -99
                        else
                            n_cells = n_cells + 1
                            lg(mn(1):mn(3),mn(2):mn(4)) = k   !! DO WE OVERWRITE AND FAVOR LARGER CELLS
                            quad_idx(mn(1):mn(3),mn(2):mn(4)) = n_cells
                            if (use_2d_flow) then
                                call set_2d_computational_lines(&
                                    n_line_u, n_line_v, k, m, n, mn, lg, lgrmin, area_mask_padded, quad_idx&
                                )
                            endif
                        endif
                    endif
                enddo
            enddo
        enddo
        deallocate(area_mask_padded)
        write(*,*) '** INFO: No. active 2D computational cells: ', n_cells
        write(*,*) '** INFO: Number of 2D Surface flow lines is: ', n_line_u, n_line_v

    end subroutine find_active_2d_comp_cells


end module m_quadtree