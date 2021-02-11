module m_lines

    type :: ClassLines
        integer :: lintot
        integer :: linall
        integer :: liutot
        integer :: livtot
        integer :: l2dtot
        integer :: l2dall
        integer :: l2dobc
        integer :: l2grbc
        integer :: l1dtot
        integer :: l1dobc
        integer, allocatable :: id(:)
        integer, allocatable :: line(:,:)
        integer, allocatable :: kcu(:)
        double precision, allocatable :: coordinates(:,:)
        double precision, allocatable :: edge_bounds(:,:)
        !double precision, allocatable :: dpumax(:)
    end type ClassLines

    integer, parameter :: TWOD_OW_LINE = 100
    integer, parameter :: TWOD_OBSTACLE_LINE = 101

    contains

    subroutine init_lines(self, nodes, model_area)

        use m_nodes

        type(ClassLines) :: self
        type(ClassNodes), intent(in) :: nodes
        integer, intent(in) :: model_area(:,:)
        
        self%liutot = 0
        self%livtot = 0
        self%l2dtot = 0
        self%l2dall = 0
        self%lintot = 0
        self%linall = 0
        self%l2dobc = 0
        self%l2grbc = 0
        self%l1dtot = 0
        self%l1dobc = 0

        !call create_2d_sw_lines2(self, nodes, quadtree_grid, model_area)
        call create_2d_sw_lines(self, nodes, model_area)
        call set_initial_2d_sw_line_attrs(self, nodes)

    end subroutine init_lines

    ! subroutine create_2d_sw_lines2(self, nodes, quadtree_grid, model_area)

    !     use MessageHandling
    !     use m_nodes
    !     use m_quadtree
    !     use m_grid_utils, only : get_lg_corners, get_pix_corners, crop_pix_coords_to_raster

    !     type(ClassLines) :: self
    !     type(ClassNodes), intent(in) :: nodes
    !     type(QuadTreeFortran), intent(in) :: quadtree_grid
    !     integer, intent(in) :: model_area(:,:)
    !     integer :: i
    !     integer :: k, m, n
    !     integer :: m0, m1, n0, n1
    !     integer :: nu, nd
    !     integer :: i0,i1,j0,j1,i2,i3,j2,j3
    !     integer :: cell_pix_corners(8)
    !     integer, allocatable :: id(:)
    !     integer, allocatable :: line(:,:)
    !     integer :: lintot
    !     integer :: nod_guess
    !     integer :: ierr

    !     nod_guess = nodes%nodtot*8
    !     !call mess(LEVEL_DEBUG, 'Allocation status of line ids: ', ierr)

    !     allocate(line(nod_guess, 2))!, stat=ierr)
    !     !call mess(LEVEL_DEBUG, 'Allocation status of line administration: ', ierr)
    !     line = 0
    !     lintot = 0
    !     !!!!! Add attributes to lines like lik, lim, lin
    !     do k=1, get_kmax(quadtree_grid)
    !         do m=1, get_mmax(quadtree_grid, k) - 1
    !             do n=1, get_nmax(quadtree_grid, k)
    !                 call get_lg_corners(k, m, n, m0, m1, n0, n1)
    !                 if (is_cell_active(quadtree_grid, (/ m0, m1 /), (/ n0 ,n1 /), k)) then
    !                     call get_pix_corners(k, m, n, get_lgrmin(quadtree_grid), i0, i1, j0, j1, i2, i3, j2, j3)
    !                     call crop_pix_coords_to_raster(model_area, i0, i1, j0, j1, i2, i3, j2, j3)
    !                     cell_pix_corners = (/ i0, i1, j0, j1, i2, i3, j2, j3 /)
    !                     ! call find_lines(quadtree_grid, model_area, cell_pix_corners, nodes%quad_coord, .TRUE., num_lines, line)

    !                     if (m0>1) then
    !                         if(get_lg(quadtree_grid, m0-1,n0) == k-1 .and. any(minval(model_area(i0-1:i0,j0:j2), 1) > 0)) then
    !                             lintot = lintot + 1
    !                             line(lintot, 1:2) = (/ nodes%quad_coord(m0-1,n0), nodes%quad_coord(m0,n0) /)
    !                         endif
    !                         if(get_lg(quadtree_grid, m0-1, n1) == k-1 .and. any(minval(model_area(i0-1:i0,j3:j1), 1) > 0)) then
    !                             lintot = lintot + 1
    !                             line(lintot, 1:2) = (/ nodes%quad_coord(m0-1,n1), nodes%quad_coord(m0,n1) /)
    !                         endif
    !                     endif

    !                     !if (m < get_mmax(quadtree_grid, k)) then
    !                         if (get_lg(quadtree_grid, m1+1,n0) == k.and.any(minval(model_area(i1:i1+1,j0:j1), 1) > 0)) then
    !                             lintot = lintot + 1
    !                             line(lintot, 1:2) = (/ nodes%quad_coord(m0,n0), nodes%quad_coord(m1+1,n0) /)
    !                         endif

    !                         if (get_lg(quadtree_grid, m1+1,n0) == k-1 .and. any(minval(model_area(i1:i1+1,j0:j2), 1) > 0)) then
    !                             lintot = lintot + 1
    !                             line(lintot, 1:2) = (/ nodes%quad_coord(m0,n0), nodes%quad_coord(m1+1,n0) /)
    !                         endif
    !                         if (get_lg(quadtree_grid, m1+1,n1) == k-1 .and. any(minval(model_area(i1:i1+1,j3:j1), 1) > 0)) then
    !                             lintot = lintot + 1
    !                             line(lintot, 1:2) = (/ nodes%quad_coord(m0,n0), nodes%quad_coord(m1+1,n1) /)
    !                         endif
    !                     !endif
    !                 endif
    !             enddo
    !         enddo
    !     enddo
        
    !     do k=1, get_kmax(quadtree_grid)
    !         do n=1, get_nmax(quadtree_grid, k) - 1
    !             do m=1, get_mmax(quadtree_grid, k)
    !                 call get_lg_corners(k, m, n, m0, m1, n0, n1)
    !                 if (is_cell_active(quadtree_grid, (/ m0, m1 /), (/ n0 ,n1 /), k)) then
    !                     call get_pix_corners(k, m, n, get_lgrmin(quadtree_grid), i0, i1, j0, j1, i2, i3, j2, j3)
    !                     call crop_pix_coords_to_raster(model_area, i0, i1, j0, j1, i2, i3, j2, j3)

    !                     if (n0>1) then
    !                         if(get_lg(quadtree_grid, m0, n0-1) == k-1 .and. any(minval(model_area(i0:i2,max(1,j0-1):j0), 2) > 0)) then
    !                             lintot = lintot + 1
    !                             line(lintot, 1:2) = (/ nodes%quad_coord(m0,n0-1), nodes%quad_coord(m0,n0) /)
    !                         endif
    !                         if(get_lg(quadtree_grid, m1, n0-1) == k-1 .and. any(minval(model_area(i3:i1,max(1,j0-1):j0), 2) > 0)) then
    !                             lintot = lintot + 1
    !                             line(lintot, 1:2) = (/ nodes%quad_coord(m1,n0-1), nodes%quad_coord(m0,n0) /)
    !                         endif
    !                     endif

    !                     !if (n < get_nmax(quadtree_grid, k)) then
    !                         if (get_lg(quadtree_grid, m0, n1+1) == k .and. any(minval(model_area(i0:i1,j1:j1+1), 2) > 0)) then
    !                             lintot = lintot + 1
    !                             line(lintot, 1:2) = (/ nodes%quad_coord(m0,n0), nodes%quad_coord(m0,n1+1) /)
    !                         endif
    !                         if (get_lg(quadtree_grid, m0, n1+1) == k-1 .and. any(minval(model_area(i0:i2,j1:j1+1), 2) > 0)) then
    !                             lintot = lintot + 1
    !                             line(lintot, 1:2) = (/ nodes%quad_coord(m0,n0), nodes%quad_coord(m0,n1+1) /)
    !                         endif
    !                         if (get_lg(quadtree_grid, m1, n1+1) == k-1 .and. any(minval(model_area(i3:i1,j1:j1+1), 2) > 0)) then
    !                             lintot = lintot + 1
    !                             line(lintot, 1:2) = (/ nodes%quad_coord(m0,n0), nodes%quad_coord(m1,n1+1) /)
    !                         endif
    !                     !endif
    !                 endif
    !             enddo
    !         enddo
    !     enddo

    !     call mess(LEVEL_INFO, 'Number of 2D Surface flow lines is: ', lintot)

    !     allocate(self%id(lintot)) 
    !     allocate(self%line(lintot, 2))
    !     self%id = (/ (i, i = 1, lintot) /)
    !     self%line = line(1:lintot, :)
    !     self%lintot = lintot
    !     deallocate(line)

    ! end subroutine create_2d_sw_lines2

    subroutine create_2d_sw_lines(self, nodes, model_area)

        use MessageHandling
        use m_nodes
        use m_grid_utils, only : get_lg_corners, get_pix_corners, crop_pix_coords_to_raster

        type(ClassLines) :: self
        type(ClassNodes), intent(in) :: nodes
        integer, intent(in) :: model_area(:,:)
        integer :: i
        integer :: nod
        integer :: nodtot
        integer :: k, m, n
        integer :: m0, m1, n0, n1
        integer :: i0,i1,j0,j1,i2,i3,j2,j3
        integer, allocatable :: line_u(:,:), line_v(:,:)
        integer :: lintot_u, lintot_v
        integer :: neighbour
        integer :: ierr
        integer, allocatable :: nodk(:)

        nodtot = get_nodtot(nodes)
        !call mess(LEVEL_DEBUG, 'Allocation status of line ids: ', ierr)

        allocate(line_u(nodtot*4, 2))!, stat=ierr)
        allocate(line_v(nodtot*4, 2))!, stat=ierr)
        !call mess(LEVEL_DEBUG, 'Allocation status of line administration: ', ierr)
        line_u = 0
        line_v = 0
        lintot_u = 0
        lintot_v = 0

        nodk = get_nodk(nodes, 1, nodtot)
        !! U direction
        do nod=1,nodtot
            k = nodk(nod)
            m = get_nodm(nodes, nod)
            n = get_nodn(nodes, nod)
            call get_lg_corners(k, m, n, m0, m1, n0, n1)
            call get_pix_corners(k, m, n, get_lgrmin(nodes), i0, i1, j0, j1, i2, i3, j2, j3)
            call crop_pix_coords_to_raster(model_area, i0, i1, j0, j1, i2, i3, j2, j3)

            !!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!! U - direction !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!
            neighbour = get_node_from_quad(nodes, m1+1,n0) 
            if (nodk(neighbour) == k.and.any(minval(model_area(i1:i1+1,j0:j1), 1) > 0)) then
                lintot_u = lintot_u + 1
                line_u(lintot_u, 1:2) = (/ nod, neighbour /)
            else
                neighbour = get_node_from_quad(nodes, m1+1,n0)
                if (nodk(neighbour) == k-1 .and. any(minval(model_area(i1:i1+1,j0:j2), 1) > 0)) then
                    lintot_u = lintot_u + 1
                    line_u(lintot_u, 1:2) = (/ nod, neighbour /)
                endif
                neighbour = get_node_from_quad(nodes, m1+1,n1)
                if (nodk(neighbour) == k-1 .and. any(minval(model_area(i1:i1+1,j3:j1), 1) > 0)) then
                    lintot_u = lintot_u + 1
                    line_u(lintot_u, 1:2) = (/ nod, neighbour /)
                endif
            endif

            neighbour = get_node_from_quad(nodes, m0-1,n0)
            if(nodk(neighbour) == k-1 .and. any(minval(model_area(i0-1:i0,j0:j2), 1) > 0)) then
                lintot_u = lintot_u + 1
                line_u(lintot_u, 1:2) = (/ neighbour, nod /)
            endif
            neighbour = get_node_from_quad(nodes, m0-1, n1)
            if(nodk(neighbour) == k-1 .and. any(minval(model_area(i0-1:i0,j3:j1), 1) > 0)) then
                lintot_u = lintot_u + 1
                line_u(lintot_u, 1:2) = (/ neighbour, nod /)
            endif      
            !!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!! U - direction !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!


            !!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!! V - direction !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!
            neighbour = get_node_from_quad(nodes, m0, n1+1)
            if (nodk(neighbour) == k .and. any(minval(model_area(i0:i1,j1:j1+1), 2) > 0)) then
                lintot_v = lintot_v + 1
                line_v(lintot_v, 1:2) = (/ nod, neighbour /)
            else
                neighbour = get_node_from_quad(nodes, m0, n1+1)
                if (nodk(neighbour) == k-1 .and. any(minval(model_area(i0:i2,j1:j1+1), 2) > 0)) then
                    lintot_v = lintot_v + 1
                    line_v(lintot_v, 1:2) = (/ nod, neighbour /)
                endif
                neighbour = get_node_from_quad(nodes, m1, n1+1)
                if (nodk(neighbour) == k-1 .and. any(minval(model_area(i3:i1,j1:j1+1), 2) > 0)) then
                    lintot_v = lintot_v + 1
                    line_v(lintot_v, 1:2) = (/ nod, neighbour /)
                endif
            endif
            neighbour = get_node_from_quad(nodes, m0, n0-1)
            if(nodk(neighbour) == k-1 .and. any(minval(model_area(i0:i2,max(1,j0-1):j0), 2) > 0)) then
                lintot_v = lintot_v + 1
                line_v(lintot_v, 1:2) = (/ neighbour, nod /)
            endif
            neighbour = get_node_from_quad(nodes, m1, n0-1)
            if(nodk(neighbour) == k-1 .and. any(minval(model_area(i3:i1,max(1,j0-1):j0), 2) > 0)) then
                lintot_v = lintot_v + 1
                line_v(lintot_v, 1:2) = (/ neighbour, nod /)
            endif
            !!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!! V - direction !!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!
        enddo

        call mess(LEVEL_INFO, 'Number of 2D Surface flow lines is: ', lintot_u+lintot_v)

        allocate(self%line(lintot_u+lintot_v, 2))
        self%line(1:lintot_u, :) = line_u(1:lintot_u, :)
        self%line(lintot_u+1:lintot_u+lintot_v, :) = line_v(1:lintot_v, :)
        self%liutot = lintot_u
        self%livtot = lintot_v
        self%l2dtot = lintot_u + lintot_v
        self%l2dall = lintot_u + lintot_v
        self%lintot = lintot_u + lintot_v
        self%linall = lintot_u + lintot_v
        deallocate(line_u, line_v)

    end subroutine create_2d_sw_lines

    subroutine set_initial_2d_sw_line_attrs(self, nodes)

        use MessageHandling
        use m_nodes, only : ClassNodes, get_nodtot, get_node_coords, get_node_bnds
        use parameters, only : NODATA

        type(ClassLines) :: self
        type(ClassNodes), intent(in) :: nodes
        double precision, allocatable :: cell_bnds(:,:)
        double precision, allocatable :: cell_coords(:,:)
        integer :: l, nodtot
        integer :: nd, nu

        allocate(self%id(self%linall))
        allocate(self%kcu(self%linall))
        allocate(self%coordinates(self%linall, 2))
        allocate(self%edge_bounds(self%linall, 4))
        self%id = 0
        self%kcu = 0
        self%coordinates = NODATA
        self%edge_bounds = NODATA

        nodtot = get_nodtot(nodes)
        cell_coords = get_node_coords(nodes, 1, nodtot)
        cell_bnds = get_node_bnds(nodes, 1, nodtot)

        do l=1,self%lintot
            self%id(l) = l
            self%kcu(l) = TWOD_OW_LINE
            nd = self%line(l,1)
            nu = self%line(l,2)
            self%coordinates(l,:) = 0.5d0 * (cell_coords(nd,:) + cell_coords(nu,:))
            self%edge_bounds(l,:) = (/ cell_bnds(nd, 3), cell_bnds(nd, 2), cell_bnds(nd, 3), cell_bnds(nd, 4) /)
        enddo

    end subroutine set_initial_2d_sw_line_attrs

    subroutine finalize_lines(self)

        use MessageHandling

        type(ClassLines) :: self

        
        deallocate(self%line)
        if(allocated(self%id)) deallocate(self%id)
        if(allocated(self%kcu)) deallocate(self%kcu)
        if(allocated(self%coordinates)) deallocate(self%coordinates)
        if(allocated(self%edge_bounds)) deallocate(self%edge_bounds)

        call mess(LEVEL_DEBUG, 'Cleanup Lines object')

    end subroutine finalize_lines


end module m_lines