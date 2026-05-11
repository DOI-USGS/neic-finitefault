module retrieve_surf_gf


   use constants, only : nny, ndis, nt, inptd
   implicit none
   integer, private :: nblock2, nx, nz, nfmax
   integer :: npt_bank, dt_bank
   real, private :: dep_max, dep_min, dep_step, dist_max, dist_min, d_step
   real, private :: grid_depth(50), grid_dist(1601)
   complex, private :: green_bank(inptd, 10, 1601, 50)


contains


   subroutine get_surf_gf_data(gf_file, gf_bank)
   implicit none
   character(len=300), intent(in) :: gf_file
   character(len=300), intent(out) :: gf_bank
   integer :: int2, int3, int4
   open(1, file=gf_file)
   read(1, *) npt_bank, dt_bank, int2, int3, int4, nfmax, nblock2
   read(1, *) dep_min, dep_step, nz
   read(1, *) dist_min, d_step, nx
   read(1, '(a)')gf_bank
   write(*,*)'GF bank name'
   write(*, *) gf_bank
   close(1)
   dep_max = dep_min + (nz-1)*dep_step
   dist_max = dist_min + (nx-1)*d_step
   end subroutine get_surf_gf_data

   
   subroutine get_surf_gf(gf_bank, d_min, d_max, z_min, z_max)
   !!
   !!  Load to memory GF bank, from a specified location. 
   !!
   implicit none
   character(len=300) :: gf_bank
   real :: d_min, d_max, z_min, z_max
   integer :: nx_b, nx_e, nz_b, nz_e, ixx, izz
   integer :: index_rec
   real :: fault_bounds(2, 2)
   real :: dt_c, t0, dep, dis0
   integer :: iz, iz0, k, k0, npt_c, ntc, n_com, nptf
   
   call grid_properties()
   
   fault_bounds(1, 1) = d_min
   fault_bounds(2, 1) = d_max
   fault_bounds(1, 2) = z_min
   fault_bounds(2, 2) = z_max
   call get_subgrid_bounds(fault_bounds, nx_b, nx_e, nz_b, nz_e)
   
   open(50, file=gf_bank, status='old', access='direct', recl=nblock2) ! open green functions bank
   
   index_rec = 0
   write(*,*) nz_b, nz_e, nx_b, nx_e
   do iz = nz_b, nz_e
      do k = nx_b, nx_e
         index_rec = (iz - 1) * nx + k
         ixx = k - nx_b + 1
         izz = iz - nz_b + 1
         read(50, rec = index_rec) iz0, k0, dis0, t0, dep, dt_c, npt_c, nptf, &
         & ((green_bank(ntc, n_com, ixx, izz), ntc = 1, nptf), n_com = 1, 10)
!
! anticipate potential errors
!                 
         if (abs(dis0 - grid_dist(k)) .gt. 1e-4) then
            write(*, *)"Error: distance to source should be ", grid_dist(k), " but real distance is ", dis0
            write(*, *)"Check size of gf bank"
            write(*, *)"Check size of 'green_bank' variable"
            write(*, *)"Check given range of distances and depths agrees with size of gf bank"
            write(*, *)"You may need to recompute GF bank"
            stop
         end if
         if (abs(dep - grid_depth(iz)) .gt. 1e-4) then
            write(*, *)"Error: depth of source should be ", grid_depth(iz), " but real depth is ", dep
            write(*, *)"Check size of gf bank"
            write(*, *)"Check size of 'green_bank' variable"
            write(*, *)"Check given range of distances and depths agrees with size of gf bank"
            write(*, *)"You may need to recompute GF bank"
            stop
         end if
      end do
   end do
   end subroutine get_surf_gf 


   function interp_gf(distance, depth, d_min, d_max, zu_min, zu_max) result(green_out)
   !!
   !!  Get GF at a given location and depth, as the average of the 4 closest
   !!  GF in the distance-depth grid. 
   !!
   implicit none
   real(8) :: distance
   real :: depth, zu_min, zu_max
   real :: fault_bounds(2, 2)
   complex :: green_out(inptd, 10)

   real :: d_min, d_max, d_min2, d_max2, z_max, z_min
   complex :: green_up_left, green_down_left, green_up_right, green_down_right
   real :: ratio_x, ratio_z, dep_use, subgrid_dist(1601), subgrid_depth(20)
   integer :: k, k_down, k_left, k_up, k_right, n, ncom
   integer :: nx_b, nx_e, nz_b, nz_e

   fault_bounds(1, 1) = d_min
   fault_bounds(2, 1) = d_max
   fault_bounds(1, 2) = zu_min
   fault_bounds(2, 2) = zu_max
   call get_subgrid_bounds(fault_bounds, nx_b, nx_e, nz_b, nz_e)
   call get_subgrid(fault_bounds, subgrid_dist, subgrid_depth)
   dep_use = depth
   green_out(:, :) = 0.0

   d_min2 = subgrid_dist(1)
   d_max2 = subgrid_dist(nx_e - nx_b + 1)
   z_min = subgrid_depth(1)
   z_max = subgrid_depth(nz_e - nz_b + 1)
   if (distance .gt. d_max2 .or. distance .lt. d_min2) then
      write(*,*)"ooh The input distance is over the x_boundary"
      write(*,*)"green function is set to zero"
      write(*,*) distance, d_max, d_min
      return
   end if
   dep_use = max(dep_use, z_min)
   dep_use = min(dep_use, z_max)
   if (depth .gt. z_max .or. depth .lt. z_min) then
      write(*, *)"The input depth is outside vertical boundary"
      write(*, *)"use the z_max or zmin instead", depth, z_max, z_min
   end if

   k_left = 1
   do k = 1, nx_e - nx_b + 1
      if (distance .lt. subgrid_dist(k)) then
         k_left = k - 1
         exit
      end if
   end do
   k_right = k_left + 1
        
   k_up = 1
   do k = 1, nz_e - nz_b + 1
      if (dep_use .lt. subgrid_depth(k)) then
         k_up = k - 1
         exit
      end if
   end do
   k_up = min(k_up, nz_e - nz_b)
   k_down = k_up + 1
!   write(*,*) distance, k_left, k_right
!   write(*,*) depth, k_up, k_down
        
   ratio_x = (distance - subgrid_dist(k_left)) / d_step
   ratio_z = (dep_use - subgrid_depth(k_up)) / dep_step

   do ncom = 1, 10
      do n = 1, inptd
         green_up_left = green_bank(n, ncom, k_left, k_up)
         green_up_right = green_bank(n, ncom, k_right, k_up)
         green_down_left = green_bank(n, ncom, k_left, k_down)
         green_down_right = green_bank(n, ncom, k_right, k_down)
         green_out(n, ncom) = (1 - ratio_z) * (1 - ratio_x) * green_up_left &
            & + ratio_z * (1 - ratio_x) * green_down_left &
            & + (1 - ratio_z) * ratio_x * green_up_right &
            & + ratio_x * ratio_z * green_down_right
      end do
   end do
   end function interp_gf


   subroutine grid_properties()
   implicit none
   integer :: k

   do k = 1, nx
      grid_dist(k) = dist_min + (k - 1) * d_step
   end do
   do k = 1, nz
      grid_depth(k) = dep_min + (k - 1) * dep_step
   end do
   end subroutine grid_properties


   subroutine get_subgrid_bounds(fault_bounds, nx_b, nx_e, nz_b, nz_e)
   implicit none
   real, intent(in) :: fault_bounds(2, 2)
   integer, intent(out) :: nx_b, nx_e, nz_b, nz_e
   real :: d_min, d_max, zu_min, zu_max

   d_min = fault_bounds(1, 1)
   d_max = fault_bounds(2, 1)
   zu_min = fault_bounds(1, 2)
   zu_max = fault_bounds(2, 2)

   if (dist_max .lt. d_max) then
      d_max = dist_max
      write(*,*)"required maximum distance is large the Green"
      write(*,*)"Function range, use the dist_max instead"
   end if

   nx_b = int((d_min - dist_min) / d_step)
   if (nx_b .le. 0) nx_b = 1
   nx_e = int((d_max - dist_min) / d_step + 0.5) + 2
   if (nx_e .gt. nx) nx_e = nx
   nz_b = int((zu_min - dep_min) / dep_step)
   if (nz_b .le. 0) nz_b = 1
   nz_e = int((zu_max - dep_min) / dep_step + 0.5) + 2
   if (nz_e .gt. nz) nz_e = nz

   end subroutine get_subgrid_bounds


   subroutine get_subgrid(fault_bounds, subgrid_dist, subgrid_depth)
   implicit none
   real, intent(in) :: fault_bounds(2, 2)
   real, intent(out) :: subgrid_depth(20), subgrid_dist(1601)
   integer :: k, nx_b, nx_e, nz_b, nz_e

   call get_subgrid_bounds(fault_bounds, nx_b, nx_e, nz_b, nz_e)

   do k = nz_b, nz_e
      subgrid_depth(k - nz_b + 1) = grid_depth(k)
   end do

   do k = nx_b, nx_e
      subgrid_dist(k - nx_b + 1) = grid_dist(k)
   end do
   end subroutine get_subgrid


end module retrieve_surf_gf
