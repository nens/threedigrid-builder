module m_quadtree

    use parameters, only : NODATA
    use iso_c_binding
    use iso_fortran_env, only : int16

    implicit none

    integer, parameter :: LINESTRING = 1
    integer, parameter :: POLY = 3
    integer, parameter :: UP = 1
    integer, parameter :: DOWN = 0

    contains

    subroutine set_refinement(refine_id, refine_geom, n0, n1, refine_level,&
        refine_type, bbox, mmax, nmax, dx, j0, lg, i0, i1) bind(c, name="f_set_refinement")
        
        use m_grid_utils, only : get_cell_geom, find_cell_intersects,&
                                 geom_in_polygon,&
                                 feature_in_bbox

        integer(kind=c_int), intent(in) :: n0
        integer(kind=c_int), intent(in) :: n1
        integer(kind=c_int), intent(in) :: j0
        integer(kind=c_int), intent(in) :: i0
        integer(kind=c_int), intent(in) :: i1
        real(kind=c_double), intent(in) :: refine_geom(n0,n1) ! 2D array (x,y) of the refinement geometry
        integer(kind=c_int), intent(in) :: refine_id ! 2D ID of the refinement geometry  
        integer(kind=c_int), intent(in) :: refine_level ! The level of k to be set by the given geometry.
        integer(kind=c_int), intent(in) :: refine_type ! Type of geometry passed in the refine_geom variable (LINESTRING or POLYGON)
        real(kind=c_double), intent(in) :: bbox(4) ! Bounding box of total Quadtree
        integer(kind=c_int), intent(in) :: mmax(j0) ! X Dimension of each refinement level
        integer(kind=c_int), intent(in) :: nmax(j0) ! Y Dimension of each refinement level
        real(kind=c_double), intent(in) :: dx(j0) ! Cell size of each refinement level
        integer(kind=c_int), intent(inout) :: lg(i0,i1) ! Array with all refinement levels.
        integer :: status
        integer :: m, n
        integer :: mnmin(2), mnmax(2)
        double precision :: cell_geom(5,2)
        double precision :: refine_geom_bbox(4)
        logical :: cross
        
        status = 0
        if(.not.refine_type==LINESTRING.and..not.refine_type==POLY) then
            write(*,*) '** WARNING: Refinement type not known. Skipping: ', refine_id, refine_type
            status = 0
            return
        endif

        refine_geom_bbox = (/ minval(refine_geom(:,1)), minval(refine_geom(:,2)),&
                                maxval(refine_geom(:,1)), maxval(refine_geom(:,2)) /)
        if (feature_in_bbox(refine_geom_bbox, bbox)) then
            mnmin = convert_to_grid_crd(bbox(1:2), dx(1), refine_geom_bbox(1:2), mmax, nmax, DOWN)
            mnmax = convert_to_grid_crd(bbox(1:2), dx(1), refine_geom_bbox(3:4), mmax, nmax, UP)
        else
            write(*,*) '** INFO: Refinement outside area_mask. ID and type: ', refine_id, refine_type
            status = 0
            return
        endif

        write(*,*) '** INFO: Start applying refinement with refinement level and type: ', refine_id, refine_level, refine_type
        do n=mnmin(2), mnmax(2)
            do m=mnmin(1), mnmax(1)
                cross=.FALSE.
                cell_geom = get_cell_geom(bbox(1), bbox(2), m, n, dx(1))
                if (minval(refine_geom(:, 1)) > minval(cell_geom(:, 1)) .and.&
                    maxval(refine_geom(:, 1)) < maxval(cell_geom(:, 1)) .and.&
                    minval(refine_geom(:, 2)) > minval(cell_geom(:, 2)).and.&
                    maxval(refine_geom(:, 1)) < maxval(cell_geom(:, 2))) then
                    cross = .TRUE.
                elseif (minval(cell_geom(:,1))>=minval(refine_geom(:,1)) - dx(1).and.&
                    maxval(cell_geom(:,1))<=maxval(refine_geom(:,1)) + dx(1).and.&
                    minval(cell_geom(:,2))>=minval(refine_geom(:,2)) - dx(1).and.&
                    maxval(cell_geom(:,2))<=maxval(refine_geom(:,2)) + dx(1)) then
                    if (refine_type==LINESTRING) then
                        cross = find_cell_intersects(refine_geom, cell_geom)
                    elseif(refine_type==POLY) then
                        cross = geom_in_polygon(refine_geom, cell_geom)
                    endif
                endif

                if (cross) then
                    lg(m,n) = min(lg(m,n), refine_level)
                    status = 1
                else
                    cycle
                endif
            enddo
        enddo

        if (status == 0) then
            write(*,*) '** WARNING: Unsuccessfully applied refinement geometry with id: ', refine_id, ' and type: ', refine_type
        endif

    end subroutine set_refinement

    subroutine make_quadtree(kmax, mmax, nmax, lgrmin, area_mask, lg, quad_idx,&
        n0, n1, i0, i1, n_cells, n_line_u, n_line_v) bind(c, name="make_quadtree")
    !!! Entry point for creating quadtree by setting lg array and then finding number of active cells and lines.

        integer(kind=c_int), intent(in) :: n0
        integer(kind=c_int), intent(in) :: n1
        integer(kind=c_int), intent(in) :: i0
        integer(kind=c_int), intent(in) :: i1
        integer(kind=c_int), intent(in) :: kmax ! Maximum refinement levels
        integer(kind=c_int), intent(in) :: mmax(kmax) ! X Dimension of each refinement level
        integer(kind=c_int), intent(in) :: nmax(kmax) ! Y Dimension of each refinement level
        integer(kind=c_int), intent(in) :: lgrmin ! Number of pixels in cell of smallest refinement level
        integer(kind=c_int16_t), intent(in) :: area_mask(n0,n1) ! Array with active pixels of model.
        integer(kind=c_int), intent(inout) :: lg(i0,i1) ! Array with all refinement levels.
        integer(kind=c_int), intent(inout) :: quad_idx(i0,i1) ! Array with idx of cell at lg refinement locations
        integer(kind=c_int), intent(inout) :: n_cells ! counter for active cells
        integer(kind=c_int), intent(inout) :: n_line_u ! counter for active u lines
        integer(kind=c_int), intent(inout) :: n_line_v ! counter for active v lines
        integer :: k
        integer :: m, n
        
        write(*,*) '** INFO: Start making quadtree.'
        do m=1, mmax(kmax)
            do n=1, nmax(kmax)
                call divide(kmax, m, n, lg)
            enddo
        enddo
        call balance_quadtree(kmax, mmax, nmax, lg)
        call find_active_2d_comp_cells(kmax, mmax, nmax, lgrmin, lg, area_mask, quad_idx, n_cells, n_line_u, n_line_v)
        write(*,*) '** INFO: Done making quadtree.'

    end subroutine make_quadtree

    recursive subroutine divide(k, m, n, lg) !ip, jp, 
    !!! Recursive subroutine to set correct refinement levels on lg refinement array.
        use m_grid_utils, only : get_lg_corners

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

    subroutine find_active_2d_comp_cells(kmax, mmax, nmax, lgrmin, lg, area_mask, quad_idx, n_cells, n_line_u, n_line_v)
    !!! Counting active cells and lines based on area_mask of active pixels.
        use m_grid_utils, only : get_lg_corners, get_pix_corners, crop_pix_coords_to_raster, pad_area_mask
        use m_cells, only : set_2d_computational_lines

        integer, intent(in) :: kmax
        integer, intent(in) :: mmax(:)
        integer, intent(in) :: nmax(:)
        integer, intent(in) :: lgrmin
        integer, intent(inout) :: lg(:,:)
        integer(kind=int16), intent(in) :: area_mask(:,:)
        integer, intent(inout) :: quad_idx(:,:)
        integer, intent(inout) :: n_line_u
        integer, intent(inout) :: n_line_v
        integer(kind=int16), allocatable:: area_mask_padded(:, :)
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
                            call set_2d_computational_lines(n_line_u, n_line_v, k, m, n, mn, lg, lgrmin, area_mask_padded, quad_idx)
                        endif
                    endif
                enddo
            enddo
        enddo
        deallocate(area_mask_padded)
        write(*,*) '** INFO: No. active 2D computational cells: ', n_cells
        write(*,*) '** INFO: Number of 2D Surface flow lines is: ', n_line_u, n_line_v

    end subroutine find_active_2d_comp_cells
    
    function convert_to_grid_crd(origin, dx, xy, mmax, nmax, round) result (mn)
    !!! Create pixel indexes (or grid coordinates) for lg array.
        double precision, intent(in) :: xy(2)
        double precision, intent(in) :: origin(2)
        double precision, intent(in) :: dx
        integer, intent(in) :: mmax(:)
        integer, intent(in) :: nmax(:)
        integer, intent(in) :: round
        integer :: mn(2)

        if (round == DOWN) then
            mn(1) = max(1, int(floor((xy(1) - origin(1) - dx) / dx)))
            mn(2) = max(1, int(floor((xy(2) - origin(2) - dx) / dx)))
        elseif (round == UP) then
            mn(1) = min(mmax(1), int(ceiling((xy(1) - origin(1)) / dx)))
            mn(2) = min(nmax(1), int(ceiling((xy(2) - origin(2)) / dx)))
        else
            write(*,*) '** ERROR: Rounding option not known: ', round
        endif

    end function convert_to_grid_crd


end module m_quadtree