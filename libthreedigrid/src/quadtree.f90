module m_quadtree

    use parameters, only : NODATA
    use iso_c_binding

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
        real(kind=c_double), intent(in) :: refine_geom(n0,n1)
        integer(kind=c_int), intent(in) :: refine_id
        integer(kind=c_int), intent(in) :: refine_level
        integer(kind=c_int), intent(in) :: refine_type
        real(kind=c_double), intent(in) :: bbox(4)
        integer(kind=c_int), intent(in) :: mmax(j0)
        integer(kind=c_int), intent(in) :: nmax(j0)
        real(kind=c_double), intent(in) :: dx(j0)
        integer(kind=c_int), intent(inout) :: lg(i0,i1)
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
            write(*,*) '** INFO: Refinement outside model_area. ID and type: ', refine_id, refine_type
            status = 0
            return
        endif

        write(*,*) '** INFO: Start applying refinement with refinement level and type: ', refine_id, refine_level, refine_type
        cross=.FALSE.
        do n=mnmin(2), mnmax(2)
            do m=mnmin(1), mnmax(1)
                cell_geom = get_cell_geom(bbox(1), bbox(2), m, n, dx(1))
                if (minval(cell_geom(:,1))>=minval(refine_geom(:,1)) - dx(1).and.&
                    maxval(cell_geom(:,1))<=maxval(refine_geom(:,1)) + dx(1).and.&
                    minval(cell_geom(:,2))>=minval(refine_geom(:,2)) - dx(1).and.&
                    maxval(cell_geom(:,2))<=maxval(refine_geom(:,2)) + dx(1)) then
                    if (refine_type==LINESTRING) then
                        cross = find_cell_intersects(refine_geom, cell_geom)  !!TODO: Check linestrings that fall within smallest cell of quadtree.
                    elseif(refine_type==POLY) then
                        cross = geom_in_polygon(refine_geom, cell_geom)
                    endif

                    if (cross) then
                        lg(m,n) = min(lg(m,n), refine_level)
                        status = 1
                    endif
                else
                    cycle
                endif
            enddo
        enddo

        if (status == 0) then
            write(*,*) '** WARNING: Unsuccessfully applied refinement geometry with id: ', refine_id, ' and type: ', refine_type
        endif

    end subroutine set_refinement

    subroutine make_quadtree(kmax, mmax, nmax, lgrmin, model_area, lg, quad_idx,&
        m0, n0, n1, i0, i1, num_active_nodes, cnt_line_u, cnt_line_v) bind(c, name="make_quadtree")
        integer(kind=c_int), intent(in) :: m0
        integer(kind=c_int), intent(in) :: n0
        integer(kind=c_int), intent(in) :: n1
        integer(kind=c_int), intent(in) :: i0
        integer(kind=c_int), intent(in) :: i1
        integer(kind=c_int), intent(in) :: kmax
        integer(kind=c_int), intent(in) :: mmax(m0)
        integer(kind=c_int), intent(in) :: nmax(m0)
        integer(kind=c_int), intent(in) :: lgrmin
        integer(kind=c_int), intent(in) :: model_area(n0,n1)
        integer(kind=c_int), intent(inout) :: lg(i0,i1)
        integer(kind=c_int), intent(inout) :: quad_idx(i0,i1)
        integer(kind=c_int), intent(inout) :: num_active_nodes
        integer(kind=c_int), intent(inout) :: cnt_line_u
        integer(kind=c_int), intent(inout) :: cnt_line_v
        integer :: k
        integer :: m, n
        
        write(*,*) '** INFO: Start making quadtree.'
        do m=1, mmax(kmax)
            do n=1, nmax(kmax)
                call divide(kmax, m, n, lg)
            enddo
        enddo
        call balance_quadtree(kmax, mmax, nmax, lg)
        call find_active_2d_comp_cells(kmax, mmax, nmax, lgrmin, lg, model_area, quad_idx, num_active_nodes, cnt_line_u, cnt_line_v)
        write(*,*) '** INFO: Done making quadtree.'

    end subroutine make_quadtree

    recursive subroutine divide(k, m, n, lg) !ip, jp, 

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

    subroutine find_active_2d_comp_cells(kmax, mmax, nmax, lgrmin, lg, model_area, quad_idx, num_active_nodes, cnt_line_u, cnt_line_v)

        use m_grid_utils, only : get_lg_corners, get_pix_corners, get_cell_bbox

        integer, intent(in) :: kmax
        integer, intent(in) :: mmax(:)
        integer, intent(in) :: nmax(:)
        integer, intent(in) :: lgrmin
        integer, intent(inout) :: lg(:,:)
        integer, intent(in) :: model_area(:,:)
        integer, intent(inout) :: quad_idx(:,:)
        integer, intent(inout) :: cnt_line_u
        integer, intent(inout) :: cnt_line_v
        integer :: k
        integer :: m,n
        integer :: mn(4)
        integer :: i0, i1, j0, j1, i2, i3, j2, j3
        integer :: num_active_nodes
        

        num_active_nodes = 0
        cnt_line_u = 0
        cnt_line_v = 0
        quad_idx = 0
        do k=kmax,1,-1
            do m=1,mmax(k)
                do n=1,nmax(k)
                    call get_pix_corners(k, m, n, lgrmin, i0, i1, j0, j1, i2, i3, j2, j3)
                    mn = get_lg_corners(k, m, n)
                    i1 = min(i1, size(model_area, 1))
                    j1 = min(j1, size(model_area, 2))
                    if (all(model_area(i0:i1, j0:j1) == 0)) then
                        lg(mn(1):mn(3),mn(2):mn(4)) = -99
                    else
                        if (any(lg(mn(1):mn(3),mn(2):mn(4)) == k)) then !! TODO: CHECK OF MODEL AREA CHECK IS NECESSARY???
                            num_active_nodes = num_active_nodes + 1
                            lg(mn(1):mn(3),mn(2):mn(4)) = k   !! DO WE OVERWRITE AND FAVOR LARGER CELLS
                            quad_idx(mn(1):mn(3),mn(2):mn(4)) = num_active_nodes
                            call count_2d_computational_lines(k, mn, lg, model_area, quad_idx, i0, i1, j0, j1, i2, i3, j2, j3, cnt_line_u, cnt_line_v)
                        endif
                    endif
                enddo
            enddo
        enddo
        write(*,*) '** INFO: No. active 2D computational cells: ', num_active_nodes
        write(*,*) '** INFO: Number of 2D Surface flow lines is: ', cnt_line_u, cnt_line_v

    end subroutine find_active_2d_comp_cells
    

    subroutine count_2d_computational_lines(k, mn, lg, model_area, quad_idx, i0, i1, j0, j1, i2, i3, j2, j3, cnt_line_u, cnt_line_v)

        use m_grid_utils, only : crop_pix_coords_to_raster

        integer, intent(in) :: k
        integer, intent(in) :: mn(4)
        integer, intent(in) :: lg(:,:)
        integer, intent(in) :: model_area(:,:)
        integer, intent(in) :: quad_idx(:,:)
        integer, intent(inout) :: i0,i1,j0,j1,i2, i3, j2, j3
        integer, intent(inout) :: cnt_line_u
        integer, intent(inout) :: cnt_line_v
        integer :: nod
        integer :: neighbour
        integer :: k_neighbour

        k_neighbour = 0
        !! U direction
        call crop_pix_coords_to_raster(model_area, i0, i1, j0, j1, i2, i3, j2, j3)

        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! U - direction !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!
        if (mn(3) < size(lg, 1)) then
            k_neighbour = lg(mn(3)+1, mn(2))
            if (k_neighbour == k.and.any(minval(model_area(i1:i1+1,j0:j1), 1) > 0)) then
                cnt_line_u = cnt_line_u + 1
            else
                if (k_neighbour == k-1 .and. any(minval(model_area(i1:i1+1,j0:j2), 1) > 0)) then
                    cnt_line_u = cnt_line_u + 1
                endif
                k_neighbour = lg(mn(3)+1,mn(4))
                if (k_neighbour == k-1 .and. any(minval(model_area(i1:i1+1,j3:j1), 1) > 0)) then
                    cnt_line_u = cnt_line_u + 1
                endif
            endif
        endif

        if (mn(1) > 1) then
            k_neighbour = lg(mn(1)-1, mn(2))
            if(k_neighbour == k-1 .and. any(minval(model_area(i0-1:i0,j0:j2), 1) > 0)) then
                cnt_line_u = cnt_line_u + 1
            endif
            k_neighbour = lg(mn(1)-1, mn(4))
            if(k_neighbour == k-1 .and. any(minval(model_area(i0-1:i0,j3:j1), 1) > 0)) then
                cnt_line_u = cnt_line_u + 1
            endif
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! U - direction !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!


        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! V - direction !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!
        if (mn(4) < size(lg, 2)) then
            k_neighbour = lg(mn(1), mn(4)+1)
            if (k_neighbour == k .and. any(minval(model_area(i0:i1,j1:j1+1), 2) > 0)) then
                cnt_line_v = cnt_line_v + 1
            else
                if (k_neighbour == k-1 .and. any(minval(model_area(i0:i2,j1:j1+1), 2) > 0)) then
                    cnt_line_v = cnt_line_v + 1
                endif
                k_neighbour = lg(mn(3), mn(4)+1)
                if (k_neighbour == k-1 .and. any(minval(model_area(i3:i1,j1:j1+1), 2) > 0)) then
                    cnt_line_v = cnt_line_v + 1
                endif
            endif
        endif
            
        if (mn(2) > 1) then
            k_neighbour = lg(mn(1), mn(2)-1)
            if(k_neighbour == k-1 .and. any(minval(model_area(i0:i2,max(1,j0-1):j0), 2) > 0)) then
                cnt_line_v = cnt_line_v + 1
            endif
            k_neighbour = lg(mn(3), mn(2)-1)
            if(k_neighbour == k-1 .and. any(minval(model_area(i3:i1,max(1,j0-1):j0), 2) > 0)) then
                cnt_line_v = cnt_line_v + 1
            endif
        endif
        !!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!! V - direction !!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine count_2d_computational_lines

    function convert_to_grid_crd(origin, dx, xy, mmax, nmax, round) result (mn)

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