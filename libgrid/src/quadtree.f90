module m_quadtree

    use MessageHandling
    use parameters, only : NODATA

    implicit none

    type :: QuadTreeFortran
        private
        double precision :: x0p
        double precision :: y0p
        double precision :: bbox(4)
        integer :: lgrmin
        integer :: kmax
        integer :: number_active_nodes
        integer, allocatable :: mmax(:)
        integer, allocatable :: nmax(:)
        double precision, allocatable :: dx(:)
        integer, allocatable :: lg(:,:)
    contains
        procedure :: init_quadtree
        procedure :: set_refinement
        procedure :: make_quadtree
        procedure :: balance_quadtree
        procedure :: set_active_2d_comp_cells
        procedure :: active_node_count
        procedure :: is_cell_active
        procedure :: get_lgrmin
        procedure :: get_kmax
        procedure :: get_mmax
        procedure :: get_nmax
        procedure :: get_dx
        procedure :: get_lg
        procedure :: get_origin
        procedure :: get_extent
        procedure :: convert_to_grid_crd
        procedure :: finalize_quadtree
    end type QuadTreeFortran

    integer, parameter :: LINESTRING = 1
    integer, parameter :: POLY = 2
    integer, parameter :: UP = 1
    integer, parameter :: DOWN = 0

    contains

    subroutine init_quadtree(self, x0p, y0p, dxp, imax, jmax, lgrmin, kmax)

        class(QuadTreeFortran) :: self
        double precision, intent(in) :: x0p
        double precision, intent(in) :: y0p
        double precision, intent(in) :: dxp
        integer, intent(in) :: imax
        integer, intent(in) :: jmax
        integer, intent(in) :: lgrmin
        integer, intent(in) :: kmax
        double precision :: cell_geom(5,2)
        integer :: k
        logical :: active_cell
        integer :: max_grid_x_pix
        integer :: max_grid_y_pix

        self%number_active_nodes = 0
        self%lgrmin = lgrmin
        self%kmax = kmax
        allocate(self%mmax(kmax), self%nmax(kmax), self%dx(kmax))
        self%mmax = 0
        self%nmax = 0
        self%dx = 0.0d0

        max_grid_x_pix = (lgrmin*2**(kmax-1)*(imax/(lgrmin*2**(kmax-1))))+lgrmin*2**(kmax-1)
        max_grid_y_pix = (lgrmin*2**(kmax-1)*(jmax/(lgrmin*2**(kmax-1))))+lgrmin*2**(kmax-1)
        do k=1,kmax
            self%mmax(k)=max_grid_x_pix/int((lgrmin*2**(k-1)))
            self%nmax(k)=max_grid_y_pix/int((lgrmin*2**(k-1)))
            self%dx(k) = lgrmin*2**(k-1) * dxp
            write(msgbuf,'(''Dimensions of grid level '' (i8)'' are mmax ''(i8)'' and nmax '' (i8))') k, self%mmax(k), self%nmax(k)
            call dbg_flush()
            write(msgbuf,'(''Cell size of grid level ''(i8)'' is ''(f12.4))') k, self%dx(k)
            call dbg_flush()
        enddo

        self%x0p = x0p
        self%y0p = y0p
        self%bbox = (/ x0p, y0p, x0p + max_grid_x_pix * dxp,  y0p + max_grid_y_pix * dxp /)
        allocate(self%lg(self%mmax(1), self%nmax(1)))
        self%lg = self%kmax
        
    end subroutine init_quadtree

    subroutine set_refinement(self, refine_id, refine_geom, refine_level, refine_type, status)
        
        use MessageHandling
        use m_grid_utils, only : get_cell_geom, find_cell_intersects,&
                                 get_lg_corners, geom_in_polygon,&
                                 feature_in_bbox

        class(QuadTreeFortran) :: self
        double precision, intent(in) :: refine_geom(:,:)
        integer, intent(in) :: refine_id
        integer, intent(in) :: refine_level
        integer, intent(in) :: refine_type
        integer, intent(out) :: status
        integer :: m, n
        integer :: mnmin(2), mnmax(2)
        double precision :: cell_geom(5,2)
        double precision :: refine_geom_bbox(4)
        logical :: cross
        
        status = 0
        if(.not.refine_type==LINESTRING.and..not.refine_type==POLY) then
            call mess(LEVEL_WARN, 'Refinement type not known. Skipping: ', refine_id, refine_type)
            status = 0
            return
        endif

        refine_geom_bbox = (/ minval(refine_geom(:,1)), minval(refine_geom(:,2)),&
                                maxval(refine_geom(:,1)), maxval(refine_geom(:,2)) /)
        if (feature_in_bbox(refine_geom_bbox, self%bbox)) then
            mnmin = convert_to_grid_crd(self, refine_geom_bbox(1:2), DOWN)
            mnmax = convert_to_grid_crd(self, refine_geom_bbox(3:4), UP)
        else
            call mess(LEVEL_INFO, 'Refinement outside model_area. ID and type: ', refine_id, refine_type)
            status = 0
            return
        endif

        call mess(LEVEL_INFO, 'Start applying refinement with refinement level and type: ', refine_id, refine_level, refine_type)
        cross=.FALSE.
        do n=mnmin(2), mnmax(2)
            do m=mnmin(1), mnmax(1)
                cell_geom = get_cell_geom(self%x0p, self%y0p, m, n, self%dx(1))
                if (minval(cell_geom(:,1))>=minval(refine_geom(:,1))-self%dx(1).and.&
                    maxval(cell_geom(:,1))<=maxval(refine_geom(:,1))+self%dx(1).and.&
                    minval(cell_geom(:,2))>=minval(refine_geom(:,2))-self%dx(1).and.&
                    maxval(cell_geom(:,2))<=maxval(refine_geom(:,2))+self%dx(1)) then
                    if (refine_type==LINESTRING) then
                        cross = find_cell_intersects(refine_geom, cell_geom)  !!TODO: Check linestrings that fall within smallest cell of quadtree.
                    elseif(refine_type==POLY) then
                        cross = geom_in_polygon(refine_geom, cell_geom)
                    endif

                    if (cross) then
                        self%lg(m,n) = min(self%lg(m,n), refine_level)
                        status = 1
                    endif
                else
                    cycle
                endif
            enddo
        enddo

        if (status == 0) then
            write(msgbuf, '(A,i8,A,i1)') 'Unsuccessfully applied refinement geometry with id: ', refine_id, ' and type: ', refine_type
            call warn_flush()
        endif

    end subroutine set_refinement

    subroutine make_quadtree(self)
        
        class(QuadTreeFortran) :: self
        integer :: k
        integer :: m, n
        
        
        do m=1, self%mmax(self%kmax)
            do n=1, self%nmax(self%kmax)
                call divide(self%kmax, m, n, self%lg)
            enddo
        enddo

        call balance_quadtree(self)

    end subroutine make_quadtree

    recursive subroutine divide(k, m, n, lg) !ip, jp, 

        use m_grid_utils, only : get_lg_corners

        integer, intent(in) :: k
        integer, intent(in) :: m
        integer, intent(in) :: n
        integer, allocatable, intent(inout) :: lg(:,:)
        integer :: k1
        integer :: m0, n0, m1, n1
        integer :: m_n, n_n

        call get_lg_corners(k, m, n, m0, m1, n0, n1)
        if(any(lg(m0:m1,n0:n1) < k)) then
            k1 = k-1
            m_n=2*m-1
            n_n=2*n-1
            lg(m0:m1,n0:n1) = min(lg(m0:m1,n0:n1), k1)
            call divide(k1, m_n, n_n, lg)
            call divide(k1, m_n+1, n_n, lg)
            call divide(k1, m_n, n_n+1, lg)
            call divide(k1, m_n+1, n_n+1, lg)
        endif

    end subroutine divide

    subroutine balance_quadtree(self)

        class(QuadTreeFortran) :: self
        integer :: k, k1
        integer :: m, n, m1, n1
        integer :: i, j

        do k=1,self%kmax-1
            do j=1, self%nmax(1)!, self%lgrmin
                do i=1,self%mmax(1)-1!, self%lgrmin
                    if(self%lg(i,j) == k .and. self%lg(i+1,j)>k+1) then
                        self%lg(i+1,j)=k+1
                    endif
                    if(self%lg(i,j) > k+1 .and. self%lg(i+1,j)==k) then
                        self%lg(i,j)=k+1
                    endif
                enddo
            enddo
            do j=1,self%nmax(1)-1
                do i=1, self%mmax(1)
                    if (self%lg(i,j)==k.and.self%lg(i,j+1)>k+1) then
                        self%lg(i,j+1)=k+1
                    endif
                    if (self%lg(i,j)>k+1.and.self%lg(i,j+1)==k) then
                        self%lg(i,j)=k+1
                    endif
                enddo
            enddo

            
            do n=1, self%nmax(self%kmax)
                do m=1, self%mmax(self%kmax)
                    call divide(self%kmax, m, n, self%lg)
                enddo
            enddo
        enddo

    end subroutine balance_quadtree


    subroutine set_active_2d_comp_cells(self, model_area)

        use MessageHandling
        use m_grid_utils, only : get_lg_corners, get_pix_corners, get_cell_bbox

        class(QuadTreeFortran) :: self
        integer, intent(in) :: model_area(:,:)
        integer :: k
        integer :: m,n
        integer :: m0, m1, n0, n1
        integer :: i0, i1, j0, j1

        self%number_active_nodes=0
        do k=1,self%kmax
            do m=1,self%mmax(k)
                do n=1,self%nmax(k)
                    call get_pix_corners(k, m, n, self%lgrmin, i0, i1, j0, j1)
                    call get_lg_corners(k, m, n, m0, m1, n0, n1)
                    i1 = min(i1, size(model_area, 1))
                    j1 = min(j1, size(model_area, 2))
                    if(all(model_area(i0:i1, j0:j1) == 0)) then
                        self%lg(m0:m1,n0:n1) = -99
                    else
                        if (any(self%lg(m0:m1,n0:n1) == k)) then !! TODO: CHECK OF MODEL AREA CHECK IS NECESSARY???
                            self%number_active_nodes = self%number_active_nodes+1
                            self%lg(m0:m1,n0:n1) = k   !! DO WE OVERWRITE AND FAVOR LARGER CELLS
                        endif
                    endif
                enddo
            enddo
        enddo
        call mess(LEVEL_INFO, 'No. active 2D computational cells: ', self%number_active_nodes)

    end subroutine set_active_2d_comp_cells

    function active_node_count(self) result(active_nodes)

        class(QuadTreeFortran) :: self
        integer :: active_nodes

        active_nodes = self%number_active_nodes

    end function active_node_count

    function is_cell_active(self, m, n, k) result(cell_active)

        class(QuadTreeFortran) :: self
        integer, intent(in) :: m(2)
        integer, intent(in) :: n(2)
        integer, intent(in) :: k
        logical :: cell_active

        cell_active = .FALSE.
        if (m(1)>0.and.m(1)<=m(2).and.m(2)<=size(self%lg,1).and.&
            n(1)>0.and.n(1)<=n(2).and.n(2)<=size(self%lg,2)) then
            if (all(self%lg(m(1):m(2),n(1):n(2)) == k)) then
                cell_active = .TRUE.
            else
                cell_active = .FALSE.
            endif
        else
            write(msgbuf, '(A,4i4)') 'Slice indices for active cell is out of bounds: ', m(1), m(2), n(1), n(2)
            call warn_flush()
        endif

    end function is_cell_active

    function get_lgrmin(self) result(lgrmin)
        
        class(QuadTreeFortran) :: self
        integer :: lgrmin

        lgrmin = self%lgrmin
    
    end function get_lgrmin

    function get_kmax(self) result(kmax)
        
        class(QuadTreeFortran) :: self
        integer :: kmax

        kmax = self%kmax
    
    end function get_kmax

    function get_mmax(self, k) result(mmax)
        
        class(QuadTreeFortran) :: self
        integer, intent(in) :: k
        integer :: mmax

        mmax = self%mmax(k)
    
    end function get_mmax

    function get_nmax(self, k) result(nmax)
        
        class(QuadTreeFortran) :: self
        integer, intent(in) :: k
        integer :: nmax

        nmax = self%nmax(k)
    
    end function get_nmax

    function get_dx(self, k) result(dx)
        
        class(QuadTreeFortran) :: self
        integer, intent(in) :: k
        double precision :: dx

        dx = self%dx(k)
    
    end function get_dx

    function get_lg(self, m0, m1, n0, n1) result(k)

        use MessageHandling

        class(QuadTreeFortran) :: self
        integer, intent(in) :: m0, m1
        integer, intent(in) :: n0, n1
        integer, allocatable :: k(:,:)

        if (m0>0.and.m0<=m1.and.m1<=size(self%lg,1).and.&
                n0>0.and.n0<=n1.and.n1<=size(self%lg,2)) then
            
            allocate(k(m0:m1,n0:n1))
            k = self%lg(m0:m1,n0:n1)
        else
            !write(msgbuf, '(A,4i4)') 'Slice indices out of bounds: ', m0, m1, n0, n1
            !call warn_flush()
            allocate(k(1,1))
            k = 0
        endif

    end function get_lg

    function get_origin(self) result(origin)

        class(QuadTreeFortran) :: self
        double precision :: origin(2)

        origin = (/ self%x0p, self%y0p /)

    end function get_origin

    function get_extent(self) result(extent)

        class(QuadTreeFortran) :: self
        double precision :: extent(4)

        extent = self%bbox

    end function get_extent

    function convert_to_grid_crd(self, xy, round) result (mn)

        use MessageHandling

        class(QuadTreeFortran) :: self
        double precision, intent(in) :: xy(2)
        integer, intent(in) :: round
        integer :: mn(2)

        if (round == DOWN) then
            mn(1) = max(1, int(floor((xy(1) - self%x0p - self%dx(1)) / self%dx(1))))
            mn(2) = max(1, int(floor((xy(2) - self%y0p - self%dx(1)) / self%dx(1))))
        elseif (round == UP) then
            mn(1) = min(self%mmax(1), int(ceiling((xy(1) - self%x0p) / self%dx(1))))
            mn(2) = min(self%nmax(1), int(ceiling((xy(2) - self%y0p) / self%dx(1))))
        else
            call mess(LEVEL_ERROR, 'Rounding option not known: ', round)
        endif

    end function convert_to_grid_crd


    subroutine finalize_quadtree(self)

        class(QuadTreeFortran) :: self

        deallocate(self%mmax)
        deallocate(self%nmax)
        deallocate(self%dx)
        if(allocated(self%lg)) deallocate(self%lg)

        call mess(LEVEL_DEBUG, 'Cleanup Quadtree object')
    
    end subroutine finalize_quadtree


end module m_quadtree