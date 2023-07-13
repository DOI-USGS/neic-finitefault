module retrieve_gf


   use constants, only : nny, ndis, nt
   implicit none
   integer, private :: block_gg, nx, nz
   integer :: lnpt_gfs
   real :: dt_gfs, t_cor
   real, private :: dep_max, dep_min, dep_step, dist_max, dist_min, d_step
   real, private :: grid_depth(nny), grid_dist(ndis)
   real, private, allocatable :: green_bank(:, :, :, :)


contains


   subroutine get_gf_data(gf_file, vel_model, gf_bank)
!
!  Args:
!  gf_file: file with properties of GF bank
!  gf_bank: location of GF bank
!  vel_model: file with velocity model
!
   implicit none
   character(len=100), intent(in) :: gf_file
   character(len=100), intent(out) :: vel_model, gf_bank
   open(1, file=gf_file, status='old')
   read(1, *)vel_model
   read(1, *)dep_max, dep_min, dep_step
   read(1, *)dist_max, dist_min, d_step
   read(1, *)lnpt_gfs, dt_gfs, block_gg, t_cor
   read(1, '(a)')gf_bank
   write(*,*)'GF bank name'
   write(*, *)gf_bank
   close(1)
   end subroutine get_gf_data

   
   subroutine get_gf(gf_bank, z_min, z_max)
!
!  Args:
!  gf_bank: location of GF bank
!  z_min: minimum depth of point sources
!  z_max: maximum depth of point sources
!
   implicit none
   character(len=100) :: gf_bank
   real :: z_min, z_max
   integer :: npt, nx_b, nx_e, nz_b, nz_e, ixx, izz
   integer :: index_rec
   real :: fault_bounds(2, 2)
   real :: dt_c, t0, dep, dis0
   integer :: iz, iz0, k, k0, npt_c, ntc, n_com
   allocate(green_bank(nt, 8, 1000, 310))
   call grid_properties()
   
   npt = 2 ** lnpt_gfs         ! careful!
  
   if (z_max .gt. dep_max) then
      write(0, *)"Error: maximum depth of GF bank is ", dep_max, " but max. depth of source is ", z_max
      write(0, *)"May need to recompute GF bank"
      stop
   endif
 
   fault_bounds(1, 1) = dist_min
   fault_bounds(2, 1) = dist_max
   fault_bounds(1, 2) = z_min
   fault_bounds(2, 2) = z_max
   call get_subgrid_bounds(fault_bounds, nx_b, nx_e, nz_b, nz_e)
   
   open(50, file=gf_bank, status='old', access='direct', recl=block_gg) ! open green functions bank
   
   index_rec = 0
   write(*,*)nz_b, nz_e, nx_b, nx_e
   do iz = nz_b, nz_e
      do k = nx_b, nx_e
         index_rec = (iz - 1) * nx + k
         ixx = k - nx_b + 1
         izz = iz - nz_b + 1
         read(50, rec = index_rec)iz0, k0, dis0, t0, dep, dt_c, npt_c, &
         & ((green_bank(ntc, n_com, ixx, izz), ntc = 1, npt), n_com = 1, 8)
!
! anticipate potential errors
!                 
         if(abs(dt_c - dt_gfs) .gt. 1e-4)then
            write(0, *)"Error: dt of GF bank file is ", dt_gfs, " but dt of data is ", dt_c
            write(0, *)"May need to recompute GF bank"
            stop
         endif
         if(abs(dis0 - grid_dist(k)) .gt. 1e-4)then
            write(0, *)"Error: distance to source should be ", grid_dist(k), " but real distance is ", dis0
            write(0, *)"Check size of gf bank"
            write(0, *)"Check size of 'green_bank' variable"
            write(0, *)"Check given range of distances and depths agrees with size of gf bank"
            write(0, *)"You may need to recompute GF bank"
            stop
         endif
         if(abs(dep - grid_depth(iz)) .gt. 1e-4)then
            write(0, *)"Error: depth of source should be ", grid_depth(iz), " but real depth is ", dis0
            write(0, *)"Check size of gf bank"
            write(0, *)"Check size of 'green_bank' variable"
            write(0, *)"Check given range of distances and depths agrees with size of gf bank"
            write(0, *)"You may need to recompute GF bank"
            stop
         endif
      enddo
   enddo
   close(50)
   end subroutine get_gf 


   function interp_gf(distance, depth, zu_min, zu_max) result(green_out)
!
!  Args:
!  distance: distance between source and station
!  depth: source depth
!  zu_min: minimum depth to get GF from GF bank
!  zu_max: maximum depth to get GF from Gf bank
!
   implicit none
   real*8 :: distance
   real :: depth, zu_min, zu_max
   real :: green_out(nt, 8), fault_bounds(2, 2)

   real :: d_min, d_max, z_max, z_min, subgrid_dist(ndis), subgrid_depth(nny)
   real :: green_up_left, green_down_left, green_up_right, green_down_right
   real :: ratio_x, ratio_z, dep_use
   integer :: k, k_down, k_left, k_up, k_right, n, ncom
   integer :: npt, nx_b, nx_e, nz_b, nz_e

   fault_bounds(1, 1) = dist_min
   fault_bounds(2, 1) = dist_max
   fault_bounds(1, 2) = zu_min
   fault_bounds(2, 2) = zu_max
   call get_subgrid_bounds(fault_bounds, nx_b, nx_e, nz_b, nz_e)
   call get_subgrid(fault_bounds, subgrid_dist, subgrid_depth)
   dep_use = depth
   green_out(:, :) = 0.0

   d_min = subgrid_dist(1)
   d_max = subgrid_dist(nx_e - nx_b + 1)
   z_min = subgrid_depth(1)
   z_max = subgrid_depth(nz_e - nz_b + 1)
   if(distance .gt. d_max .or. distance .lt. d_min)then
      write(*,*)"ooh The input distance is over the x_boundary"
      write(*,*)"green function is set to zero"
      write(*,*)distance,d_max,d_min
      return
   endif
   dep_use = max(dep_use, z_min)
   dep_use = min(dep_use, z_max)
   if(depth.gt.z_max .or. depth.lt.z_min)then
      write(*, *)"The input depth is outside vertical boundary"
      write(*, *)"use the z_max or zmin instead",depth,z_max,z_min
   endif

   k_left = 1
   do k = 1, nx_e - nx_b + 1
      if(distance .lt. subgrid_dist(k))then
         k_left = k - 1
         exit
      endif
   enddo
   k_right = k_left + 1
        
   k_up = 1
   do k = 1, nz_e - nz_b + 1
      if(dep_use .lt. subgrid_depth(k))then
         k_up = k - 1
         exit
      endif
   enddo
   k_up = min(k_up, nz_e - nz_b)
   k_down = k_up + 1
        
   ratio_x = (distance - subgrid_dist(k_left)) / d_step
   ratio_z = (dep_use - subgrid_depth(k_up)) / dep_step

   npt = 2 ** lnpt_gfs         ! careful!
   do ncom=1,8
      do n=1,npt
         green_up_left = green_bank(n, ncom, k_left, k_up)
         green_up_right = green_bank(n, ncom, k_right, k_up)
         green_down_left = green_bank(n, ncom, k_left, k_down)
         green_down_right = green_bank(n, ncom, k_right, k_down)
         green_out(n, ncom) = (1 - ratio_z) * (1 - ratio_x) * green_up_left &
            & + ratio_z * (1 - ratio_x) * green_down_left &
            & + (1 - ratio_z) * ratio_x * green_up_right &
            & + ratio_x * ratio_z * green_down_right
      enddo
   enddo
   end function interp_gf


   subroutine grid_properties()
   implicit none
   integer :: k

   nx = int((dist_max - dist_min) / d_step) + 1
   nz = int((dep_max - dep_min) / dep_step) + 1

   do k = 1, nx
      grid_dist(k) = dist_min + (k - 1) * d_step
   enddo
   do k = 1, nz
      grid_depth(k) = dep_min + (k - 1) * dep_step
   enddo
   end subroutine grid_properties


   subroutine get_subgrid_bounds(fault_bounds, nx_b, nx_e, nz_b, nz_e)
!
!  Args:
!  fault_bounds: minimum and maximum depths and distances for the fault
!  nx_b: index of minimum source-station distance
!  nx_e: index of maximum source-station distance
!  nz_b: index of minimum depth
!  nz_e: index of maximum depth
!
   implicit none
   real, intent(in) :: fault_bounds(2, 2)
   integer, intent(out) :: nx_b, nx_e, nz_b, nz_e
   real :: d_min, d_max, zu_min, zu_max

   d_min = fault_bounds(1, 1)
   d_max = fault_bounds(2, 1)
   zu_min = fault_bounds(1, 2)
   zu_max = fault_bounds(2, 2)

   if(dist_max.lt.d_max)then
      d_max=dist_max
      write(*,*)"required maximum distance is large the Green"
      write(*,*)"Function range, use the dist_max instead"
   endif

   nx_b = int((d_min - dist_min) / d_step)
   if(nx_b.le.0)nx_b = 1
   nx_e = int((d_max - dist_min) / d_step + 0.5) + 2
   if(nx_e.gt.nx)nx_e = nx
   nz_b = int((zu_min - dep_min) / dep_step)
   if(nz_b.le.0)nz_b = 1
   nz_e = int((zu_max - dep_min) / dep_step + 0.5) + 2
   if(nz_e.gt.nz)nz_e = nz

   end subroutine get_subgrid_bounds


   subroutine get_subgrid(fault_bounds, subgrid_dist, subgrid_depth)
!
!  Args:
!  fault_bounds: location of GF bank
!  subgrid_dist: range of source-station distances to be used
!  subgrid_depth: range of depths to be used
!
   implicit none
   real, intent(in) :: fault_bounds(2, 2)
   real, intent(out) :: subgrid_depth(nny), subgrid_dist(ndis)
   integer :: k, nx_b, nx_e, nz_b, nz_e

   call get_subgrid_bounds(fault_bounds, nx_b, nx_e, nz_b, nz_e)

   do k = nz_b, nz_e
      subgrid_depth(k - nz_b + 1) = grid_depth(k)
   enddo

   do k = nx_b, nx_e
      subgrid_dist(k - nx_b + 1) = grid_dist(k)
   enddo
   end subroutine get_subgrid


   subroutine deallocate_gf()
   implicit none
   deallocate(green_bank)
   end subroutine deallocate_gf


end module retrieve_gf
