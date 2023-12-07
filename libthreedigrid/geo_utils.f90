module m_grid_utils

    implicit none

    contains

    function feature_in_bbox(feature_bbox, bbox) result(in_bbox)

        double precision, intent(in) :: feature_bbox(4)
        double precision, intent(in) :: bbox(4)
        logical :: in_bbox

        if((feature_bbox(1)>=bbox(1).and.feature_bbox(3)>=bbox(1)).or.&
            (feature_bbox(1)<=bbox(3).and.feature_bbox(3)<=bbox(1)).or.&
            (feature_bbox(2)>=bbox(2).and.feature_bbox(4)>=bbox(2)).or.&
            (feature_bbox(2)<=bbox(4).and.feature_bbox(4)<=bbox(4))) then
            in_bbox = .TRUE.
        else
            in_bbox = .FALSE.
        endif

    end function feature_in_bbox


    subroutine get_pix_corners(k, m, n, min_pix, i0, i1, j0, j1, i2, i3, j2, j3)

        integer(KIND=4), intent(in) :: k
        integer, intent(in) :: m
        integer, intent(in) :: n
        integer, intent(in) :: min_pix
        integer, intent(out) :: i0
        integer, intent(out) :: i1
        integer, intent(out) :: j0
        integer, intent(out) :: j1
        integer, intent(out), optional :: i2
        integer, intent(out), optional :: i3
        integer, intent(out), optional :: j2
        integer, intent(out), optional :: j3
        integer :: num_pix

        num_pix = min_pix*2**(k-1)

        i0 = (m-1)*num_pix + 1
        i1 = m*num_pix
        j0 = (n-1)*num_pix + 1
        j1 = n*num_pix

        if (present(i2)) then
            i2 = (m-1)*num_pix + int(0.5*num_pix)
            i3 = (m-1)*num_pix + int(0.5*num_pix) + 1
            j2 = (n-1)*num_pix + int(0.5*num_pix)
            j3 = (n-1)*num_pix + int(0.5*num_pix) + 1
        endif

    end subroutine get_pix_corners

    subroutine crop_pix_coords_to_raster(raster, i0, i1, j0, j1, i2, i3, j2, j3)

        integer, intent(in) :: raster(:,:)
        integer, intent(inout) :: i0
        integer, intent(inout) :: i1
        integer, intent(inout) :: j0
        integer, intent(inout) :: j1
        integer, intent(inout), optional :: i2
        integer, intent(inout), optional :: i3
        integer, intent(inout), optional :: j2
        integer, intent(inout), optional :: j3

        i0 = min(i0, size(raster, 1))
        j0 = min(j0, size(raster, 2))
        i1 = min(i1, size(raster, 1))
        j1 = min(j1, size(raster, 2))
        i2 = min(i2, size(raster, 1))
        j2 = min(j2, size(raster, 2))
        i3 = min(i3, size(raster, 1))
        j3 = min(j3, size(raster, 2))

    end subroutine crop_pix_coords_to_raster

    function pad_area_mask(raster, i0, i1, j0, j1) result(padded_raster)

        integer*1, intent(in) :: raster(:,:)
        integer, intent(in) :: i0
        integer, intent(in) :: i1
        integer, intent(in) :: j0
        integer, intent(in) :: j1
        integer*1, allocatable :: padded_raster(:, :)
        integer :: i_size, j_size, size_raster_i, size_raster_j

        size_raster_i = size(raster, 1)
        size_raster_j = size(raster, 2)
        i_size = max(0, i1, size_raster_i)
        j_size = max(0, j1, size_raster_j)
        allocate(padded_raster(1:i_size, 1:j_size))
        padded_raster = 0

        padded_raster(1:size_raster_i, 1:size_raster_j) = raster(1:size_raster_i, 1:size_raster_j)

    end function pad_area_mask

    function get_lg_corners(k, m, n) result(mn)

        integer, intent(in) :: k
        integer, intent(in) :: m
        integer, intent(in) :: n
        integer :: mn(4)
        integer :: num_cells
        
        num_cells = 1 * 2 ** (k - 1)
        mn(1) = (m-1)*num_cells + 1
        mn(2) = (n-1)*num_cells + 1
        mn(3) = m*num_cells
        mn(4) = n*num_cells

    end function get_lg_corners

    function get_cell_bbox(xorig, yorig, m, n, dx) result(bbox)

        double precision, intent(in) :: xorig
        double precision, intent(in) :: yorig
        integer, intent(in) :: m
        integer, intent(in) :: n
        double precision, intent(in) :: dx
        double precision :: bbox(4)

        bbox(1) = xorig + m * dx - dx
        bbox(2) = yorig + n * dx - dx
        bbox(3) = xorig + m * dx
        bbox(4) = yorig + n * dx

    end function get_cell_bbox

    
    function get_cell_geom(xorig, yorig, m, n, dx) result(cell_geom)

        double precision, intent(in) :: xorig
        double precision, intent(in) :: yorig
        integer, intent(in) :: m
        integer, intent(in) :: n
        double precision, intent(in) :: dx
        double precision :: x0, y0, x1, y1
        double precision :: cell_geom(5,2)

        x0 = xorig + m * dx - dx
        y0 = yorig + n * dx - dx
        x1 = xorig + m * dx
        y1 = yorig + n * dx

        cell_geom(1,:) = (/ x0, y0 /)
        cell_geom(2,:) = (/ x0, y1 /)
        cell_geom(3,:) = (/ x1, y1 /)
        cell_geom(4,:) = (/ x1, y0 /)
        cell_geom(5,:) = (/ x0, y0 /)
    
    end function get_cell_geom

    function find_cell_intersects(geom1, geom2) result (cross)

        double precision, intent(in) :: geom1(:,:)
        double precision, intent(in) :: geom2(:,:)
        double precision :: edge_geom(4)
        logical :: cross
        integer :: i, j
        double precision :: intersection(2)

        do i=1, size(geom1, 1) - 1
            do j=1, size(geom2, 1) - 1
                cross = .FALSE.
                intersection = -9999.0d0

                edge_geom = (/ geom2(j,1), geom2(j,2), geom2(j+1, 1) , geom2(j+1, 2) /)    
                call get_line_intersection(&
                    edge_geom, (/geom1(i, 1), geom1(i, 2), geom1(i+1, 1), geom1(i+1, 2) /), cross, intersection&
                )
                if (cross.eqv..TRUE.) then
                    !write(*,*) 'intersection x, y: ', intersection
                    return
                endif
            enddo
        enddo

    end function find_cell_intersects

    subroutine find_crossing(geom1,geom2,cross,xcross,ycross)
        
        implicit none
        
        double precision, intent(in) :: geom1(4), geom2(4)
        double precision :: x0,y0,x1,y1,x2,y2,x3,y3
        logical, intent(out) :: cross
        double precision, intent(out), optional :: xcross, ycross
        double precision :: aa0,aa1,bb0,bb1,cc0,cc1,ddi,ds0,ds1,ds2
        
        !> details First it finds the corresponding connection.
        !! Second, the orientation of levee with respect to grid is determined
        !! Third, It is determined which velocity points are connected and which flooding treshold it belongs.
        !> @param[in] x0, y0 coordinates of beginning of levee element
        !> @param[in] x1, y1 coordinates of end of levee element
        !> @param[in] x2, y2 start coordinates of the second line
        !> @param[in] x3, y3 end coordinates of the second line
        x0 = geom1(1)
        y0 = geom1(2)
        x1 = geom1(3)
        y1 = geom1(4)
        x2 = geom2(1)
        y2 = geom2(2)
        x3 = geom2(3)
        y3 = geom2(4)

        cross=.FALSE.
        if ((y3.eq.y2).and.(y0.gt.y2).and.(y1.gt.y2)) return
        if ((x3.eq.x2).and.(x0.gt.x2).and.(x1.gt.x2)) return
        if ((y3.eq.y2).and.(y0.lt.y2).and.(y1.lt.y2)) return
        if ((x3.eq.x2).and.(x0.lt.x2).and.(x1.lt.x2)) return
        if ((max(y0,y1).lt.min(y2,y3)).or.(min(y0,y1).gt.max(y2,y3))) return
        if ((max(x0,x1).lt.min(x2,x3)).or.(min(x0,x1).gt.max(x2,x3))) return
        
        aa0=y1-y0
        bb0=x0-x1
        cc0=x0*y1-x1*y0
        
        aa1=y3-y2
        bb1=x2-x3
        cc1=x2*y3-x3*y2

        ddi=aa0*bb1-bb0*aa1

        xcross=(cc0*bb1-bb0*cc1)/ddi
        ycross=(aa0*cc1-cc0*aa1)/ddi
        ds0=dsqrt((x3-x2)**2+(y3-y2)**2)
        ds1=dsqrt((x3-xcross)**2+(y3-ycross)**2)
        ds2=dsqrt((xcross-x2)**2+(ycross-y2)**2)
        
        if (ds0-1.0d-8>ds1+ds2.or.abs(ddi)>1.0d-8) then ! points are not on the same side of the levee or lines are not parallel
            cross = .TRUE.
        endif
        
    end subroutine find_crossing


    subroutine find_intersection(geom1, geom2, intersects, intersection)

        implicit none

        double precision, intent(in) :: geom1(4)
        double precision, intent(in) :: geom2(4)
        logical, intent(out) :: intersects
        double precision, intent(out), optional :: intersection(2)
        double precision :: a1
        double precision :: a2
        double precision :: b1
        double precision :: b2
        !real      :: x, y       ! intersect point
        double precision :: dx1
        double precision :: dx2
        double precision :: dy1
        double precision :: dy2
       
        intersects = .FALSE.
        if ((max(geom1(2),geom1(4)).lt.min(geom2(2),geom2(4))).or.(min(geom1(2),geom1(4)).gt.max(geom2(2),geom2(4)))) then
            return
        endif
        if ((max(geom1(1),geom1(3)).lt.min(geom2(1),geom2(3))).or.(min(geom1(1),geom1(3)).gt.max(geom2(1),geom2(3)))) then
            return
        endif

        dx1 = geom1(1) - geom1(3)
        dy1 = geom1(2) - geom1(4)
        if( dx1 == 0.0d0) then    ! in case this line is of the form y = b
            a1 = 0.0d0
            b1 = geom1(2)
        else
            a1= dy1 / dx1
            b1 = geom1(2) - a1*geom1(1)
        endif

        dx2 = geom2(1) - geom2(3)
        dy2 = geom2(2) - geom2(4)
        if(dx2 == 0.0d0) then    ! in case this line is of the form y = b
            a2 = 0.0d0
            b2 = geom2(2)
        else
            a2= dy2 / dx2
            b2 = geom2(2) - a2*geom2(1)
        endif
        
        !write(*,*) 'a1', a1, 'a2', a2, 'b1', b1, 'b2', b2

        if( (a1 - a2)<1.0d-5 ) then
            intersects = .FALSE.
            write(*,*) '** DEBUG: Lines do not intersect'
            return
        else
            intersects = .TRUE.
            if(present(intersection)) then
                intersection(1) = (b2-b1) / (a1-a2)
                intersection(2) = a1 * intersection(1) + b1
            endif
            return
        endif
    
    end subroutine find_intersection

    subroutine get_line_intersection(geom1, geom2, intersect, intersection)
        
        implicit none

        double precision, intent(in) :: geom1(4)
        double precision, intent(in) :: geom2(4)
        logical, intent(out) :: intersect
        double precision, intent(out), optional :: intersection(2)
        double precision :: s1_x, s1_y, s2_x, s2_y
        double precision :: s, t, det, eps

        intersect = .FALSE.

        if ((geom2(4).eq.geom2(2)).and.(geom1(2).gt.geom2(2)).and.(geom1(4).gt.geom2(2))) return
        if ((geom2(3).eq.geom2(1)).and.(geom1(1).gt.geom2(1)).and.(geom1(3).gt.geom2(1))) return
        if ((geom2(4).eq.geom2(2)).and.(geom1(2).lt.geom2(2)).and.(geom1(4).lt.geom2(2))) return
        if ((geom2(3).eq.geom2(1)).and.(geom1(1).lt.geom2(1)).and.(geom1(3).lt.geom2(1))) return
        if ((max(geom1(2),geom1(4)).lt.min(geom2(2),geom2(4))).or.(min(geom1(2),geom1(4)).gt.max(geom2(2),geom2(4)))) then
            return
        endif
        if ((max(geom1(1),geom1(3)).lt.min(geom2(1),geom2(3))).or.(min(geom1(1),geom1(3)).gt.max(geom2(1),geom2(3)))) then
            return
        endif

        eps = 1.0d-5
        s1_x = geom1(3) - geom1(1)
        s1_y = geom1(4) - geom1(2)
        s2_x = geom2(3) - geom2(1)
        s2_y = geom2(4) - geom2(2)

        det = -s2_x * s1_y + s1_x * s2_y
        if(det<=eps.and.det>=-eps) then !! Lines are linear
            return
        endif
        s = (-s1_y * (geom1(1) - geom2(1)) + s1_x * (geom1(2) - geom2(2))) / det
        t = ( s2_x * (geom1(2) - geom2(2)) - s2_y * (geom1(1) - geom2(1))) / det

        if (s >= 0 .and. s <= 1 .and. t >= 0 .and. t <= 1) then
        !! Collision detected
            if (present(intersection)) then
                intersection(1) = geom1(1) + (t * s1_x);
                intersection(2) = geom1(2) + (t * s1_y);
            endif
            intersect = .TRUE.
            return
        endif
        
    end subroutine get_line_intersection

    function geom_in_polygon(geom1, geom2) result(cross)

        double precision, intent(in) :: geom1(:,:)
        double precision, intent(in) :: geom2(:,:)
        integer :: in_polygon
        logical :: cross
        integer :: i

        cross = .FALSE.
        do i=1,size(geom2,1)-1
            call point_in_polygon(geom2(i,1), geom2(i,2), geom1(:,1), geom1(:,2), in_polygon)
            if (in_polygon==0.or.in_polygon==1) then
                cross = .TRUE.
                !write(*,*) 'intersection x, y: ', intersection
                return
            endif
        enddo

    end function geom_in_polygon

    subroutine point_in_polygon(p_x, p_y, poly_x, poly_y, in_polygon)

        implicit none

        double precision, intent(in) :: p_x
        double precision, intent(in) :: p_y
        double precision, intent(in) :: poly_x(:)
        double precision, intent(in) :: poly_y(:)
        integer :: i, j, n_vertices
        integer, intent(out) :: in_polygon
        
        double precision :: xi, yi, xj, yj
        logical :: ix, iy, jx, jy, eor
	
    !   EXCLUSIVE OR STATEMENT FUNCTION.
        eor(ix, iy) = (ix.or.iy) .and. .not.(ix.and.iy)
	
        in_polygon = -1
        n_vertices = size(poly_x)
	
        do i= 1, n_vertices
            xi = poly_x(i) - p_x
            yi = poly_y(i) - p_y
	!CHECK WHETHER THE POINT IN QUESTION IS AT THIS VERTEX.
            if( xi==0.0d0.and.yi==0.0d0) then
                in_polygon = 0
                return
            endif
	!J IS NEXT VERTEX NUMBER OF POLYGON.
            j = 1 + mod(i,n_vertices)
            xj = poly_x(j) - p_x
            yj = poly_y(j) - p_y
	!IS THIS LINE OF 0 LENGTH ?
            if(xi==xj.and.yi==yj ) then
                cycle
            endif
            ix = xi>=0.0d0
            iy = yi>=0.0d0
            jx = xj>=0.0d0
            jy = yj>=0.0d0
	!CHECK WHETHER (PX,PY) IS ON VERTICAL SIDE OF POLYGON.
            if(xi==0.0d0.and.xj==0.0d0.and.eor(iy,jy)) then
                in_polygon = 0
                return
            endif
	!CHECK WHETHER (PX,PY) IS ON HORIZONTAL SIDE OF POLYGON.
            if(yi==0.0d0.and.yj==0.0d0.and.eor(ix,jx)) then
                in_polygon = 0
                return
            endif
	!CHECK WHETHER BOTH ENDS OF THIS SIDE ARE COMPLETELY 1) TO RIGHT
	!OF, 2) TO LEFT OF, OR 3) BELOW (PX,PY).
            if(.not.((iy.or.jy).and.eor(ix,jx))) then
                cycle
            endif
	!DOES THIS SIDE OBVIOUSLY CROSS LINE RISING VERTICALLY FROM (PX,PY)
            if(.not.(iy.and.jy.and.eor(ix,jx))) then
                if((yi*xj-xi*yj)/(xj-xi)<0.0d0) then
                    cycle
                elseif((yi*xj-xi*yj)/(xj-xi)==0.0d0) then
                    in_polygon = 0
                    return
                else
                    in_polygon = -in_polygon
                endif
            else
                in_polygon = -in_polygon
            endif
	
        enddo

    end subroutine point_in_polygon


end module m_grid_utils
