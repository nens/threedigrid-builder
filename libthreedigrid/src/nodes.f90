module m_nodes

    type :: ClassNodes
        private
        integer :: nodtot
        integer :: nodall
        integer :: n2dtot
        integer :: n2dall
        integer :: n2dobc
        integer :: n2grbc
        integer :: n1dtot
        integer :: n1dobc
        integer :: lgrmin
        integer, allocatable :: id(:)
        integer, allocatable :: type(:)
        integer, allocatable :: nodk(:)
        integer, allocatable :: nodm(:)
        integer, allocatable :: nodn(:)
        integer, allocatable :: quad_coord(:,:)

        double precision, allocatable :: coordinates(:,:)
        double precision, allocatable :: bounds(:,:)
        double precision, allocatable :: sumax(:)
    contains
        procedure :: init_nodes
        procedure :: init_2d_node_attributes
        procedure :: map_grid_to_nodes
        procedure :: get_node_coords
        procedure :: get_id
        procedure :: get_node_bnds
        procedure :: get_nodtot
        procedure :: get_nodk
        procedure :: get_nodm
        procedure :: get_nodn
        procedure :: get_lgrmin
        procedure, private :: get_node_from_quad_scalar
        procedure, private :: get_node_from_quad_slice
        procedure :: get_nodgrid
        procedure :: finalize_nodes
        generic   :: get_node_from_quad =>  get_node_from_quad_scalar,&
                                            get_node_from_quad_slice
    end type ClassNodes

    ! interface get_node_from_quad
        ! module procedure get_node_from_quad_scalar
        ! module procedure get_node_from_quad_slice
    ! end interface get_node_from_quad


    contains

    subroutine init_nodes(self, quadtree_grid)

        use m_quadtree, only : QuadTreeFortran, active_node_count,&
                                get_lgrmin, get_mmax, get_nmax
                                
        class(ClassNodes) :: self
        class(QuadTreeFortran), intent(in) :: quadtree_grid
        integer :: shape_lg(2)

        self%n2dobc = 0
        self%n2grbc = 0
        self%n1dtot = 0
        self%n1dobc = 0
        self%n2dtot = quadtree_grid%active_node_count()
        self%n2dall = quadtree_grid%active_node_count()
        self%nodtot = quadtree_grid%active_node_count()
        self%nodall = quadtree_grid%active_node_count()
        self%lgrmin = quadtree_grid%get_lgrmin()

        shape_lg = (/ quadtree_grid%get_mmax(1), quadtree_grid%get_nmax(1) /)
        call self%init_2d_node_attributes(shape_lg)
        call self%map_grid_to_nodes(quadtree_grid)

    end subroutine init_nodes

    
    subroutine init_2d_node_attributes(self, shape_lg)

        use parameters, only : NODATA

        class(ClassNodes) :: self
        integer, intent(in) :: shape_lg(2)

        allocate(self%id(self%nodtot))
        allocate(self%type(self%nodtot))
        allocate(self%nodk(self%nodtot))
        allocate(self%nodm(self%nodtot))
        allocate(self%nodn(self%nodtot))
        allocate(self%quad_coord(shape_lg(1),shape_lg(2)))
        allocate(self%coordinates(self%nodtot ,2))
        allocate(self%bounds(self%nodtot, 4))

        self%id = 0
        self%type = 0
        self%nodk = 0
        self%nodm = 0
        self%nodn = 0
        self%quad_coord = -99
        self%coordinates = NODATA
        self%bounds = NODATA

    end subroutine init_2d_node_attributes

    
    subroutine map_grid_to_nodes(self, quadtree_grid)

        use MessageHandling
        use m_quadtree, only : QuadTreeFortran,&
                               get_kmax, get_mmax, get_nmax,&
                               get_dx, is_cell_active, get_origin
        use m_grid_utils, only : get_lg_corners, get_cell_bbox
        use parameters, only : NODE_2D,&
                               NODE_2D_GW,&
                               NODE_1D,&
                               NODE_1D_STOR,&
                               NODE_2D_BOUND,&
                               NODE_2D_GW_BOUND,&
                               NODE_1D_BOUND,&
                               NODATA

        class(ClassNodes) :: self
        class(QuadTreeFortran), intent(in) :: quadtree_grid

        integer :: nod
        integer :: k
        integer :: m, n
        integer :: m0, m1, n0, n1
        double precision :: origin(2)
        double precision :: bbox(4)

        nod = 0
        do k=quadtree_grid%get_kmax(),1,-1
            do m=1,quadtree_grid%get_mmax(k)
                do n=1,quadtree_grid%get_nmax(k)
                    call get_lg_corners(k, m, n, m0, m1, n0, n1)
                    if (quadtree_grid%is_cell_active((/ m0, m1 /), (/ n0 ,n1 /), k)) then
                        nod = nod + 1
                        self%id(nod) = nod
                        self%type(nod) = NODE_2D
                        self%nodk(nod) = k
                        self%nodm(nod) = m
                        self%nodn(nod) = n
                        self%quad_coord(m0:m1,n0:n1) = nod                        
                        origin = quadtree_grid%get_origin()
                        bbox = get_cell_bbox(origin(1), origin(2), m, n, quadtree_grid%get_dx(k))
                        self%bounds(nod,:) = bbox
                        self%coordinates(nod, :) = (/ 0.5d0 * (bbox(1) + bbox(3)), 0.5d0 * (bbox(2) + bbox(4)) /)
                    else
                        continue
                    endif
                enddo
            enddo
        enddo
        call mess(LEVEL_INFO, 'Number of 2D nodes is: ', nod)

    end subroutine map_grid_to_nodes

    function get_node_coords(self, n_start, n_end) result(xy)
        
        use m_array_utils, only : check_bounds

        class(ClassNodes), target :: self
        integer, intent(in) :: n_start
        integer, intent(in) :: n_end
        double precision, pointer :: xy(:,:)

        if (check_bounds(self%coordinates, n_start, n_end)) then
            xy => self%coordinates(n_start:n_end, 1:2)
        else
            xy => NULL()
        endif

    end function get_node_coords

    function get_id(self, n_start, n_end) result(id)
        
        use m_array_utils, only : check_bounds

        class(ClassNodes), target :: self
        integer, intent(in) :: n_start
        integer, intent(in) :: n_end
        integer, pointer :: id(:)

        if (check_bounds(self%id, n_start, n_end)) then
            id => self%id(n_start:n_end)
        else
            id => NULL()
        endif

    end function get_id

    function get_node_bnds(self, n_start, n_end) result(bounds)

        use m_array_utils, only : check_bounds

        class(ClassNodes), target :: self
        integer, intent(in) :: n_start
        integer, intent(in) :: n_end
        double precision, pointer :: bounds(:,:)
        
        if (check_bounds(self%bounds, n_start, n_end)) then
            bounds => self%bounds(n_start:n_end, 1:4)
            write(*,*) 'SHAPE: ', shape(bounds)
        else
            bounds => NULL()
        endif

    end function get_node_bnds

    function get_nodtot(self) result(nodtot)

        class(ClassNodes) :: self
        integer :: nodtot

        nodtot = self%nodtot

    end function get_nodtot

    function get_nodk(self, n_start, n_end) result(k)

        class(ClassNodes) :: self
        integer, intent(in) :: n_start
        integer, intent(in) :: n_end
        integer :: k(1:n_end-n_start)

        k = self%nodk(n_start:n_end)

    end function get_nodk

    function get_nodm(self, nod) result(m)

        class(ClassNodes) :: self
        integer, intent(in) :: nod
        integer :: m

        m = self%nodm(nod)

    end function get_nodm

    function get_nodn(self, nod) result(n)

        class(ClassNodes) :: self
        integer, intent(in) :: nod
        integer :: n

        n = self%nodn(nod)

    end function get_nodn

    function get_lgrmin(self) result(lgrmin)

        class(ClassNodes) :: self
        integer :: lgrmin

        lgrmin = self%lgrmin

    end function get_lgrmin

    function get_node_from_quad_scalar(self, m, n) result(nod)

        class(ClassNodes) :: self
        integer, intent(in) :: m
        integer, intent(in) :: n
        integer :: nod

        if (m>0.and.m<=size(self%quad_coord,1).and.n>0.and.n<=size(self%quad_coord,2)) then
            nod = self%quad_coord(m,n)
        else
            !call mess(LEVEL_WARN, 'Slice indices out of bounds: ', m, n)
            nod = 0
        endif

    end function get_node_from_quad_scalar

    function get_node_from_quad_slice(self, m0, m1, n0, n1) result(nod)

        use MessageHandling

        class(ClassNodes) :: self
        integer, intent(in) :: m0, m1
        integer, intent(in) :: n0, n1
        integer, allocatable :: nod(:,:)

        if (m0>0.and.m0<=m1.and.m1<=size(self%quad_coord,1).and.&
                n0>0.and.n0<=n1.and.n1<=size(self%quad_coord,2)) then
            
            allocate(nod(m0:m1,n0:n1))
            nod = self%quad_coord(m0:m1,n0:n1)
        else
            !write(msgbuf, '(A,4i4)') 'Slice indices out of bounds: ', m0, m1, n0, n1
            !call warn_flush()
            allocate(nod(1,1))
            nod = 0
        endif

    end function get_node_from_quad_slice

    subroutine get_nodgrid(self, mask, nodgrid)

        use m_grid_utils, only : get_pix_corners

        class(ClassNodes) :: self
        logical, intent(in) :: mask(:,:)
        integer, intent(inout) :: nodgrid(:,:)
        integer :: nod
        integer :: ip0, ip1, jp0, jp1
        integer :: size_i, size_j
        integer :: i, j
        
        size_i = size(nodgrid, 1)
        size_j = size(nodgrid, 2)

        do nod=1,self%n2dtot
            call get_pix_corners(self%nodk(nod), self%nodm(nod), self%nodn(nod), self%lgrmin, ip0, ip1, jp0, jp1)
            ip0 = max(1, ip0)
            jp0 = max(1, jp0)

            ip1 = min(size_i, ip1)
            jp1 = min(size_j, jp1)
            if (ip0>=1.and.ip1<=+size_i.and.&
                    jp0>=1.and.jp1<=+size_j) then
                
                where(mask(ip0:ip1,jp0:jp1).eqv..TRUE.) nodgrid(ip0:ip1,jp0:jp1)=nod

            endif
        enddo

    end subroutine get_nodgrid

    subroutine finalize_nodes(self)

        use MessageHandling

        class(ClassNodes) :: self

        deallocate(self%id)
        deallocate(self%type)
        deallocate(self%nodk)
        deallocate(self%nodm)
        deallocate(self%nodn)
        deallocate(self%quad_coord)
        deallocate(self%coordinates)
        deallocate(self%bounds)

        call mess(LEVEL_DEBUG, 'Cleanup Nodes object')
        
    end subroutine finalize_nodes

end module m_nodes