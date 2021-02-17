! module m_lines

!     type :: ClassLines
!         integer :: lintot
!         integer :: linall
!         integer :: liutot
!         integer :: livtot
!         integer :: l2dtot
!         integer :: l2dall
!         integer :: l2dobc
!         integer :: l2grbc
!         integer :: l1dtot
!         integer :: l1dobc
!         integer, allocatable :: id(:)
!         integer, allocatable :: line(:,:)
!         integer, allocatable :: kcu(:)
!         double precision, allocatable :: coordinates(:,:)
!         double precision, allocatable :: edge_bounds(:,:)
!         !double precision, allocatable :: dpumax(:)
!     contains
!         procedure :: init_lines
!         procedure :: create_2d_sw_lines
!         procedure :: set_initial_2d_sw_line_attrs
!         procedure :: finalize_lines
!     end type ClassLines

!     integer, parameter :: TWOD_OW_LINE = 100
!     integer, parameter :: TWOD_OBSTACLE_LINE = 101

!     contains

!     subroutine init_lines(self, nodes, model_area)

!         use m_nodes

!         class(ClassLines) :: self
!         class(ClassNodes), intent(in) :: nodes
!         integer, intent(in) :: model_area(:,:)
        
!         self%liutot = 0
!         self%livtot = 0
!         self%l2dtot = 0
!         self%l2dall = 0
!         self%lintot = 0
!         self%linall = 0
!         self%l2dobc = 0
!         self%l2grbc = 0
!         self%l1dtot = 0
!         self%l1dobc = 0

!         !call create_2d_sw_lines2(self, nodes, quadtree_grid, model_area)
!         call self%create_2d_sw_lines(nodes, model_area)
!         call self%set_initial_2d_sw_line_attrs(nodes)

!     end subroutine init_lines

!     subroutine create_2d_sw_lines(self, nodes, model_area)

!         use MessageHandling
!         use m_nodes
!         use m_grid_utils, only : get_lg_corners, get_pix_corners, crop_pix_coords_to_raster

!         class(ClassLines) :: self
!         class(ClassNodes), intent(in) :: nodes
!         integer, intent(in) :: model_area(:,:)
!         integer :: i
!         integer :: nod
!         integer :: nodtot
!         integer :: k, m, n
!         integer :: m0, m1, n0, n1
!         integer :: i0,i1,j0,j1,i2,i3,j2,j3
!         integer, allocatable :: line_u(:,:), line_v(:,:)
!         integer :: lintot_u, lintot_v
!         integer :: neighbour
!         integer :: ierr
!         integer, allocatable :: nodk(:)

!         nodtot = nodes%get_nodtot()
!         !call mess(LEVEL_DEBUG, 'Allocation status of line ids: ', ierr)

!         allocate(line_u(nodtot*4, 2))!, stat=ierr)
!         allocate(line_v(nodtot*4, 2))!, stat=ierr)
!         !call mess(LEVEL_DEBUG, 'Allocation status of line administration: ', ierr)
!         line_u = 0
!         line_v = 0
!         lintot_u = 0
!         lintot_v = 0

!         nodk = nodes%get_nodk(1, nodtot)
!         !! U direction
!         do nod=1,nodtot
!             k = nodk(nod)
!             m = nodes%get_nodm(nod)
!             n = nodes%get_nodn(nod)
!             call get_lg_corners(k, m, n, m0, m1, n0, n1)
!             call get_pix_corners(k, m, n, nodes%get_lgrmin(), i0, i1, j0, j1, i2, i3, j2, j3)
!             call crop_pix_coords_to_raster(model_area, i0, i1, j0, j1, i2, i3, j2, j3)

!             !!!!!!!!!!!!!!!!!!!!!!!!!
!             !!!!! U - direction !!!!!
!             !!!!!!!!!!!!!!!!!!!!!!!!!
!             neighbour = nodes%get_node_from_quad(m1+1,n0) 
!             if (nodk(neighbour) == k.and.any(minval(model_area(i1:i1+1,j0:j1), 1) > 0)) then
!                 lintot_u = lintot_u + 1
!                 line_u(lintot_u, 1:2) = (/ nod, neighbour /)
!             else
!                 neighbour = nodes%get_node_from_quad(m1+1,n0)
!                 if (nodk(neighbour) == k-1 .and. any(minval(model_area(i1:i1+1,j0:j2), 1) > 0)) then
!                     lintot_u = lintot_u + 1
!                     line_u(lintot_u, 1:2) = (/ nod, neighbour /)
!                 endif
!                 neighbour = nodes%get_node_from_quad(m1+1,n1)
!                 if (nodk(neighbour) == k-1 .and. any(minval(model_area(i1:i1+1,j3:j1), 1) > 0)) then
!                     lintot_u = lintot_u + 1
!                     line_u(lintot_u, 1:2) = (/ nod, neighbour /)
!                 endif
!             endif

!             neighbour = nodes%get_node_from_quad(m0-1,n0)
!             if(nodk(neighbour) == k-1 .and. any(minval(model_area(i0-1:i0,j0:j2), 1) > 0)) then
!                 lintot_u = lintot_u + 1
!                 line_u(lintot_u, 1:2) = (/ neighbour, nod /)
!             endif
!             neighbour = nodes%get_node_from_quad(m0-1, n1)
!             if(nodk(neighbour) == k-1 .and. any(minval(model_area(i0-1:i0,j3:j1), 1) > 0)) then
!                 lintot_u = lintot_u + 1
!                 line_u(lintot_u, 1:2) = (/ neighbour, nod /)
!             endif      
!             !!!!!!!!!!!!!!!!!!!!!!!!!
!             !!!!! U - direction !!!!!
!             !!!!!!!!!!!!!!!!!!!!!!!!!


!             !!!!!!!!!!!!!!!!!!!!!!!!!
!             !!!!! V - direction !!!!!
!             !!!!!!!!!!!!!!!!!!!!!!!!!
!             neighbour = nodes%get_node_from_quad(m0, n1+1)
!             if (nodk(neighbour) == k .and. any(minval(model_area(i0:i1,j1:j1+1), 2) > 0)) then
!                 lintot_v = lintot_v + 1
!                 line_v(lintot_v, 1:2) = (/ nod, neighbour /)
!             else
!                 neighbour = nodes%get_node_from_quad(m0, n1+1)
!                 if (nodk(neighbour) == k-1 .and. any(minval(model_area(i0:i2,j1:j1+1), 2) > 0)) then
!                     lintot_v = lintot_v + 1
!                     line_v(lintot_v, 1:2) = (/ nod, neighbour /)
!                 endif
!                 neighbour = nodes%get_node_from_quad(m1, n1+1)
!                 if (nodk(neighbour) == k-1 .and. any(minval(model_area(i3:i1,j1:j1+1), 2) > 0)) then
!                     lintot_v = lintot_v + 1
!                     line_v(lintot_v, 1:2) = (/ nod, neighbour /)
!                 endif
!             endif
!             neighbour = nodes%get_node_from_quad(m0, n0-1)
!             if(nodk(neighbour) == k-1 .and. any(minval(model_area(i0:i2,max(1,j0-1):j0), 2) > 0)) then
!                 lintot_v = lintot_v + 1
!                 line_v(lintot_v, 1:2) = (/ neighbour, nod /)
!             endif
!             neighbour = nodes%get_node_from_quad(m1, n0-1)
!             if(nodk(neighbour) == k-1 .and. any(minval(model_area(i3:i1,max(1,j0-1):j0), 2) > 0)) then
!                 lintot_v = lintot_v + 1
!                 line_v(lintot_v, 1:2) = (/ neighbour, nod /)
!             endif
!             !!!!!!!!!!!!!!!!!!!!!!!!!
!             !!!!! V - direction !!!!!
!             !!!!!!!!!!!!!!!!!!!!!!!!!
!         enddo

!         call mess(LEVEL_INFO, 'Number of 2D Surface flow lines is: ', lintot_u+lintot_v)

!         allocate(self%line(lintot_u+lintot_v, 2))
!         self%line(1:lintot_u, :) = line_u(1:lintot_u, :)
!         self%line(lintot_u+1:lintot_u+lintot_v, :) = line_v(1:lintot_v, :)
!         self%liutot = lintot_u
!         self%livtot = lintot_v
!         self%l2dtot = lintot_u + lintot_v
!         self%l2dall = lintot_u + lintot_v
!         self%lintot = lintot_u + lintot_v
!         self%linall = lintot_u + lintot_v
!         deallocate(line_u, line_v)

!     end subroutine create_2d_sw_lines

!     subroutine set_initial_2d_sw_line_attrs(self, nodes)

!         use MessageHandling
!         use m_nodes, only : ClassNodes, get_nodtot, get_node_coords, get_node_bnds
!         use parameters, only : NODATA

!         class(ClassLines) :: self
!         class(ClassNodes), intent(in) :: nodes
!         double precision, pointer :: cell_bnds(:,:)
!         double precision, pointer :: cell_coords(:,:)
!         integer :: l, nodtot
!         integer :: nd, nu

!         allocate(self%id(self%linall))
!         allocate(self%kcu(self%linall))
!         allocate(self%coordinates(self%linall, 2))
!         allocate(self%edge_bounds(self%linall, 4))
!         self%id = 0
!         self%kcu = 0
!         self%coordinates = NODATA
!         self%edge_bounds = NODATA

!         nodtot = nodes%get_nodtot()
!         cell_coords => nodes%get_node_coords(1, nodtot)
!         cell_bnds => nodes%get_node_bnds(1, nodtot)

!         do l=1,self%lintot
!             self%id(l) = l
!             self%kcu(l) = TWOD_OW_LINE
!             nd = self%line(l,1)
!             nu = self%line(l,2)
!             self%coordinates(l,:) = 0.5d0 * (cell_coords(nd,:) + cell_coords(nu,:))
!             self%edge_bounds(l,:) = (/ cell_bnds(nd, 3), cell_bnds(nd, 2), cell_bnds(nd, 3), cell_bnds(nd, 4) /)
!         enddo

!     end subroutine set_initial_2d_sw_line_attrs

!     function get_lintot(self) result(lintot)

!         class(ClassLines) :: self
!         integer :: lintot

!         lintot = self%lintot

!     end function get_lintot

!     function get_id(self, l0, l1) result(id)
        
!         use m_array_utils, only : check_bounds

!         class(ClassLines), target :: self
!         integer, intent(in) :: l0
!         integer, intent(in) :: l1
!         integer, pointer :: id(:)

!         if (check_bounds(self%id, l0, l1)) then
!             id => self%id(l0:l1)
!         else
!             id => NULL()
!         endif

!     end function get_id


!     function get_line(self, l0, l1) result(line)

!         use m_array_utils, only : check_bounds
        
!         class(ClassLines), target :: self
!         integer, intent(in) :: l0
!         integer, intent(in) :: l1
!         integer, pointer :: line(:,:)
        
!         if (check_bounds(self%line, l0, l1)) then
!             line => self%line(l0:l1, 1:2)
!         else
!             line => NULL()
!         endif

!     end function get_line

!     subroutine finalize_lines(self)

!         use MessageHandling

!         class(ClassLines) :: self

        
!         deallocate(self%line)
!         if(allocated(self%id)) deallocate(self%id)
!         if(allocated(self%kcu)) deallocate(self%kcu)
!         if(allocated(self%coordinates)) deallocate(self%coordinates)
!         if(allocated(self%edge_bounds)) deallocate(self%edge_bounds)

!         call mess(LEVEL_DEBUG, 'Cleanup Lines object')

!     end subroutine finalize_lines


! end module m_lines