module model_parameters


   use constants, only : max_seg, max_stk_psources, max_dip_psources, max_stk_subfaults, &
               &     max_dip_subfaults, max_subf, max_subfaults, max_subfaults2, max_psources
   implicit none
   real :: slip0(max_subfaults), rake0(max_subfaults), rupt_time0(max_subfaults)
   real :: t_rise0(max_subfaults), t_fall0(max_subfaults)
   integer :: segments, nxs_sub(max_seg), nys_sub(max_seg), nx_p, ny_p, nxs0, nys0
   real :: dip(max_seg), strike(max_seg), delay_seg(max_seg)
   real, allocatable :: point_sources(:, :, :)
   real :: shear(max_subfaults)
   real :: c_depth, dxs, dys, v_ref, v_min, v_max, tbl, tbr
   real :: ta0, dta
   real :: time_min(max_subfaults), time_max(max_subfaults), rake_min
   real :: time_ref(max_subfaults)
   integer :: msou
   real :: beg(max_subfaults2), dp(max_subfaults2)
   integer :: np(max_subfaults2)
   integer :: subfaults, cum_subfaults(max_seg)
   integer :: event_segment(max_seg)
!
! for regularization
!
   integer :: nleft(3, max_subfaults), nright(3, max_subfaults), & 
   & nup(3, max_subfaults), ndown(3, max_subfaults)


contains


   subroutine get_faults_data()
   implicit none
   integer io_x, io_y, iys, io_v_d, k, psource, nxys(max_seg), &
   &  segment, ll, io_seg, subfault
   real :: dip_s, stk_s, shear2(max_subf, max_seg)
   allocate(point_sources(7, max_psources, max_subfaults))
!
!     Input Fault position to memory
!
   write(*,*)'Read and store fault segments data to memory...'
   open(12, file='fault&rise_time.txt', status='old')
   read(12,*) nxs0, nys0, c_depth
   read(12,*) segments, dxs, dys, nx_p, ny_p, v_min, v_max, tbl, tbr
   read(12,*) ta0, dta, msou, v_ref, io_v_d
   subfault = 0
   cum_subfaults(:) = 0
   do segment = 1, segments
      cum_subfaults(segment) = subfault
      read(12,*) io_seg, dip(segment), strike(segment)
      read(12,*) nxs_sub(segment), nys_sub(segment), delay_seg(segment)
      do ll = 1, nxs_sub(segment) * nys_sub(segment)
         subfault = subfault + 1
         read(12,*) slip0(subfault), rake0(subfault), &
         & rupt_time0(subfault), t_rise0(subfault), t_fall0(subfault)
!  magic
         slip0(subfault) = int(slip0(subfault))
         rake0(subfault) = int(rake0(subfault))
!
      end do
   end do
   close(12)
   subfaults = subfault
!
!     Input Fault model to memory
!
   open(12, file='point_sources.txt', status='old')
   subfault = 0
   do segment = 1, segments
      read(12,*) io_seg, dip_s, stk_s
      if ((abs(dip_s-dip(segment)) .gt. 1.e-2).or. &
      &  (abs(stk_s-strike(segment)) .gt. 1.e-2)) then
         write(*,*)'Fault mechanism in Fault.pos is not matched with that in Fault.das'
         write(*,*) segment
         write(*,*) dip_s, dip(segment)
         write(*,*) stk_s, strike(segment)
      end if

      do IYS = 1, NYS_sub(segment)*NXS_sub(segment)
         subfault = subfault+1
         do psource = 1, ny_p*nx_p
            read(12,*)(point_sources(k, psource, subfault), k = 1, 7)
         end do
      end do
   end do
   close(12)
!  
!  Input shear modulous Model
!
   open(12, file='shear_model.txt', status='old')
   read(12,*) io_x
   if (io_x.ne.segments) then
      write(*,*)'Amount of fault segments is different between mu file and Fault file'
      stop
   end if
   subfault = 0
   do segment = 1, segments
      read(12,*) io_seg, io_x, io_y
      nxys(segment) = io_x*io_y
      read(12,*)(shear2(k, segment), k = 1, nxys(segment))
      do k = 1, nxys(segment)
         subfault = subfault + 1
         shear(subfault) = shear2(k, segment) 
      enddo
   end do 
   close(12)
   end subroutine get_faults_data


   subroutine events_segments()
   implicit none
   integer :: integer1, i
   open(20, file='segments_events.txt', status='old')
   read(20,*)
   do i=1, segments
      read(20, *) integer1, event_segment(i)
   enddo
   close(20)
   end subroutine events_segments


   subroutine get_model_space()
   implicit none
!
!  boundary conditions? or what?
!
   integer :: nx0, ny0
   parameter(nx0 = max_stk_subfaults*2, ny0 = max_dip_subfaults*2)
   integer :: io_surf, nblock, npa, k, segment, io_right, io_left, io_up, io_down 
   real :: delt_x, delt_y, zmed_max, zmed_min, zleft_max, zleft_min
   real :: zright_max, zright_min, zup_max, zup_min, zdown_max, zdown_min
   real :: angle_max, angle_min, vel_max, vel_min
   real :: ddx1, ddx2, ddy1, ddy2
   integer :: io_seg, nmed, nleft2, nright2, nup2, ndown2, nangle, npv, nb, nsour
   integer :: nx, ny, i, j, k_s, i_ss, i_x, i_y
   real :: surface(1000, 4), xyb(nx0, ny0, 3), xr(5), u0(max_subfaults2, 4)
   open(17, file='model_space.txt', status='old')
   read(17,*) io_surf
   if (io_surf .eq. 1) then
      open(16, file='surface.constrain', status='old')
      read(16,*) nblock, io_surf
      do i = 1, nblock
         read(16,*)(surface(i, k), k = 1, 4)
      end do
      close(16)
   end if
   npa = 0
   k = 0
   do segment = 1, segments
      read(17,*) io_seg
      read(17,*) delt_x, delt_y
      read(17,*) io_left, io_right, io_up, io_down
      read(17,*) zmed_max, zmed_min, nmed
      read(17,*) zleft_max, zleft_min, nleft2
      read(17,*) zright_max, zright_min, nright2
      read(17,*) zup_max, zup_min, nup2
      read(17,*) zdown_max, zdown_min, ndown2
      read(17,*) angle_max, angle_min, nangle
      read(17,*) vel_max, vel_min, npv
      read(17,*) nb, nsour
      nx = nxs_sub(segment)
      ny = nys_sub(segment)
      rake_min = angle_min
!
!  check whether the contrains is frighting eath other
!
      if (io_right .eq. 1 .or. io_left .eq. 1) then
         ddx1 = (zmed_max-zleft_max)/(int((nx+1)/2))
         ddx2 = (zmed_max-zright_max)/(int((nx+1)/2))
         if (delt_x .lt. ddx1 .or. delt_x .lt. ddx2) then
            write(*,*)'the constrain of left or right boundary is error'
            stop
         end if
      end if
      if (io_down .eq. 1 .or. io_up .eq. 1) then
         ddy1 = (zmed_max-zup_max)/(int((ny+1)/2))
         ddy2 = (zmed_max-zdown_max)/(int((ny+1)/2))
         if (delt_y .lt. ddy1 .or. delt_y .lt. ddy2) then
            write(*,*)'the constrain of left or right boundary is error'
            stop
         end if
      end if

!  First:  Give the range of mediate part of space

      do i = 2, nxs_sub(segment)-1
         do j = 2, nys_sub(segment)-1
            xyb(i, j, 1) = zmed_min
            xyb(i, j, 2) = zmed_max
            xyb(i, j, 3) = nmed
         end do
      end do
!
!     Second: Give the value range of left and right part
!
      do j = 1, nys_sub(segment)
         if (io_left .eq. 1) then
            xyb(1, j, 1) = zleft_min
            xyb(1, j, 2) = zleft_max
            xyb(1, j, 3) = nleft2
         else
            xyb(1, j, 1) = zmed_min
            xyb(1, j, 2) = zmed_max
            xyb(1, j, 3) = nmed
         end if
         if (io_right .eq. 1) then
            xyb(nxs_sub(segment), j, 1) = zright_min
            xyb(nxs_sub(segment), j, 2) = zright_max
            xyb(nxs_sub(segment), j, 3) = nright2
         else
            xyb(nxs_sub(segment), j, 1) = zmed_min
            xyb(nxs_sub(segment), j, 2) = zmed_max
            xyb(nxs_sub(segment), j, 3) = nmed
         end if
      end do
!
!  finally: Give the value range of up and down part
!
      do i = 2, nxs_sub(segment)-1
         if (io_up .eq. 1) then
            xyb(i, 1, 1) = zup_min
            xyb(i, 1, 2) = zup_max
            xyb(i, 1, 3) = nup2
         else
            xyb(i, 1, 1) = zmed_min
            xyb(i, 1, 2) = zmed_max
            xyb(i, 1, 3) = nmed
         end if
         if (io_down .eq. 1) then
            xyb(i, nys_sub(segment), 1) = zdown_min
            xyb(i, nys_sub(segment), 2) = zdown_max
            xyb(i, nys_sub(segment), 3) = ndown2
         else
            xyb(i, nys_sub(segment), 1) = zmed_min
            xyb(i, nys_sub(segment), 2) = zmed_max
            xyb(i, nys_sub(segment), 3) = nmed
         end if
      end do
!
!  Recheck the range of mediate part
!
      do i = 2, nxs_sub(segment)-1
         do j = 2, nys_sub(segment)-1
            xr(1) = (i-1)*delt_x+zleft_max
            xr(2) = (nxs_sub(segment)-i)*delt_x+zright_max
            xr(3) = (j-1)*delt_y+zup_max
            xr(4) = (nys_sub(segment)-j)*delt_y+zdown_max
            xr(5) = xyb(i, j, 2)
            call bbsort(xr, 1, 5)
            xyb(i, j, 2) = xr(1)
         end do
      end do
!  
!  check the surface constrain
!  
      if (io_surf .eq. 1) then
         do k_s = 1, nblock
            i_ss = int(surface(k_s, 1)+0.1)
            if (i_ss .eq. segment) then
               i_y = 1
               i_x = int(surface(k_s, 2)+0.1)
               xyb(i_x, i_y, 1) = surface(k_s, 3)
               xyb(i_x, i_y, 2) = surface(k_s, 4)
            end if
         end do
      end if
!
!  Change xy_range to xyb    
!
      npa = npa+4*nxs_sub(segment)*nys_sub(segment)
      do j = 1, nys_sub(segment)
         do i = 1, nxs_sub(segment)
            k = k+1
            u0(k, 1) = xyb(i, j, 1)
            u0(k, 2) = xyb(i, j, 2)
            u0(k, 3) = xyb(i, j, 3)
            k = k+1
            u0(k, 1) = angle_min
            u0(k, 2) = angle_max
            u0(k, 3) = nangle
            k = k+1
            u0(k, 1) = vel_min
            u0(k, 2) = vel_max
            u0(k, 3) = npv
            k = k+1
            u0(k, 1) = 1
            u0(k, 2) = MSOU
            u0(k, 3) = MSOU
         end do
      end do
   end do
   close(17)

   do k = 1, npa
      beg(k) = u0(k, 1)
      np(k) = int(u0(k, 3) + 0.1)
      if (np(k) .gt. 1) then
         dp(k) = (u0(k, 2)-u0(k, 1))/(np(k)-1)
      else
         dp(k) = 0
      end if
   end do
   end subroutine get_model_space

   
   subroutine bbsort(a, mm, nn)
!
!  Args:
!  a: array of values
!  mm: 
!  nn: 
!  
   implicit none
   real :: a(:), d
   integer :: mm, nn, m, j, i
   m = nn-mm+1
   do while (m .gt. 0)
! 10      if (m .gt. 0) then
      j = m+mm-2
      m = 0
      do i = mm, j
         if (a(i) .gt. a(i+1)) then
            d = a(i)
            a(i) = a(i+1)
            a(i+1) = d
            m = i-mm+1
         end if
      end do
   end do
!        go to 10
!        end if
   end subroutine bbsort


   subroutine get_special_boundaries()
   implicit none
   integer :: nm, nn, iss, ixs, iys, ll, segment, kxy
!   real :: dd(max_subf, max_seg), aa(max_subf, max_seg)
!
!  special boundary
!
   open(12, file='special_model_space.txt', status='old')
   read(12,*) nm
   do nn = 1, nm
      read(12,*) iss, ixs, iys
      ll = 0
      do segment = 1, iss-1
         ll = ll+nxs_sub(segment)*nys_sub(segment)
      end do
      kxy = ixs+(iys-1)*nxs_sub(iss)
      ll = ll+ixs+(iys-1)*nxs_sub(iss)
      np(4*(ll-1)+1) = 2
      dp(4*(ll-1)+1) = 10
      beg(4*(11-1)+1) = 1
!      dd(kxy, iss) = beg(4*(ll-1)+1)
      np(4*(ll-1)+2) = 2
!      aa(kxy, iss) = beg(4*(ll-1)+2)
      np(4*(ll-1)+3) = 2
      np(4*(ll-1)+4) = 2
   end do
   close(12)
   end subroutine get_special_boundaries

   
   subroutine subfault_positions()
   implicit none
   integer k, nn_use
   integer n_is, subfault, subfault2
   integer ix, iy, nxx, nyy, io_change
   character(len=80) aaaa
   integer i_s
!
!  We detect adjacent subfaults and their relative location
!  
   subfault = 0
   do i_s = 1, segments
      do iy = 1, nys_sub(i_s)
         do ix = 1, nxs_sub(i_s)
            subfault = subfault + 1
            nup(1, subfault) = i_s
            nup(2, subfault) = ix
            nup(3, subfault) = iy-1
            ndown(1, subfault) = i_s
            ndown(2, subfault) = ix
            ndown(3, subfault) = iy+1
            nleft(1, subfault) = i_s
            nleft(2, subfault) = ix-1
            nleft(3, subfault) = iy
            nright(1, subfault) = i_s
            nright(2, subfault) = ix+1
            nright(3, subfault) = iy
         end do
      end do
   end do
   open(22, file='regularization_borders.txt', status='old')
   do I_s = 1, segments 
      subfault = cum_subfaults(i_s)
!      up 
      read(22,*)
      read(22,*) n_is, nyy
      subfault2 = subfault
      do ix = 1, nxs_sub(I_s)
         subfault2 = subfault2 + 1
         nup(1, subfault2) = n_is
         nup(2, subfault2) = ix
         nup(3, subfault2) = nyy
      end do
!      down
      read(22,*) n_is, nyy
      nn_use = nys_sub(i_s)
      subfault2 = subfault + (nys_sub(i_s)-1)*nxs_sub(i_s)
      do ix = 1, nxs_sub(I_s)
         subfault2 = subfault2 + 1
         ndown(1, subfault2) = n_is
         ndown(2, subfault2) = ix
         ndown(3, subfault2) = nyy
      end do
!      left
      read(22,*) n_is, nxx
      subfault2 = subfault + 1
      do iy = 1, nys_sub(I_s)
         nleft(1, subfault2) = n_is
         nleft(2, subfault2) = nxx
         nleft(3, subfault2) = iy
         subfault2 = subfault2 + nxs_sub(i_s)
      end do
!      right
      read(22,*) n_is, nxx
      nn_use = nxs_sub(i_s)
      subfault2 = subfault
      do iy = 1, nys_sub(I_s)
         subfault2 = subfault2 + nxs_sub(i_s)
         nright(1, subfault2) = n_is
         nright(2, subfault2) = nxx
         nright(3, subfault2) = iy
      end do
   end do
   close(22)   
   open(22, file='special_regularization_borders.txt', status='old')
   read(22,*) nn_use
   do k = 1, nn_use
      read(22,'(a)')aaaa
! 1  format(a)
!      write(*,'(a)')aaaa
      read(22,*) i_s, ix, iy
      read(22,*) io_change, n_is, nxx, nyy
      subfault = cum_subfaults(i_s) + nxs_sub(i_s)*iy + ix
      if (io_change .eq. 1) then
         nup(1, subfault) = n_is
         nup(2, subfault) = nxx
         nup(3, subfault) = nyy
      end if
      read(22,*) io_change, n_is, nxx, nyy
      if (io_change .eq. 1) then
         ndown(1, subfault) = n_is
         ndown(2, subfault) = nxx
         ndown(3, subfault) = nyy
      end if
      read(22,*) io_change, n_is, nxx, nyy
      if (io_change .eq. 1) then
         nleft(1, subfault) = n_is
         nleft(2, subfault) = nxx
         nleft(3, subfault) = nyy
      end if
      read(22,*) io_change, n_is, nxx, nyy
      if (io_change .eq. 1) then
         nright(1, subfault) = n_is
         nright(2, subfault) = nxx
         nright(3, subfault) = nyy
      end if
   end do
   close(22)
   end subroutine subfault_positions
   

   subroutine get_rise_time(ta00, dta0, msou0)
!
!  Args:
!  ta00: first parameter of rise time function
!  dta0: second parameter of rise time function
!  msou: time windows for rise time function
!  
   implicit none
   real :: ta00, dta0
   integer :: msou0 
   ta00 = ta0
   dta0 = dta
   msou0 = msou
   end subroutine get_rise_time

   
   subroutine get_shear(shear0)
!
!  Args:
!  shear0: shear modulous for each subfault
!  
   implicit none
   real :: shear0(:) 
   shear0(:) = shear(:)
   end subroutine get_shear

   
   subroutine get_segments(nxs_sub0, nys_sub0, dip0, strike0, delay_seg0, &
   &  segments0, subfaults0, cum_subfaults0, c_depth0)
!
!  Args:
!  nxs_sub0: amount of subfaults in strike dimension for each segment
!  nys_sub0: amount of subfaults in dip dimension for each segment
!  strike0: strike value for each segment
!  dip0: dip value for each segment
!  delay_seg0: time delay for each segment
!  segments0: amount of fault segments
!  subfaults0: total amount of subfaults
!  cum_subfaults0: cumulative amount of subfaults for each segment
!  c_depth0: hypocenter depth
!  
   implicit none
   integer :: nxs_sub0(:), nys_sub0(:), subfaults0, segments0, cum_subfaults0(:)
   real :: dip0(:), strike0(:), delay_seg0(:), c_depth0
   nxs_sub0(:) = nxs_sub(:)
   nys_sub0(:) = nys_sub(:)
   dip0(:) = dip(:)
   strike0(:) = strike(:)
   delay_seg0(:) = delay_seg(:)
   segments0 = segments
   subfaults0 = subfaults
   cum_subfaults0(:) = cum_subfaults(:)
   c_depth0 = c_depth
   end subroutine get_segments


   subroutine get_subfaults(dxs0, dys0, nx_p0, ny_p0, v_min0, v_max0, v_ref0)
!
!  Args:
!  dxs0: length along strike of each subfault
!  dys0: length along dip of each subfault
!  nx_p0: amount of point sources in strike dimension for each subfault
!  ny_p0: amount of point sources in dip dimension for each subfault
!  v_min0: mimimum rupture velocity allowed
!  v_max0: maximum rupture velocity allowed
!  v_ref0: reference rupture velocity
!  
   implicit none
   integer :: nx_p0, ny_p0
   real :: dxs0, dys0, v_min0, v_max0, v_ref0
   nx_p0 = nx_p
   ny_p0 = ny_p
   dxs0 = dxs
   dys0 = dys
   v_min0 = v_min
   v_max0 = v_max
   v_ref0 = v_ref
   end subroutine get_subfaults
   

   subroutine get_space(time_min0, time_max0, time_ref0, beg0, dp0, np0)
!
!  Args:
!  time_min0: minimum allowed time of rupture initiation at each subfault
!  time_max0: maximum allowed time of rupture initiation at each subfault
!  time_ref0: reference time of rupture initiation at each subfault
!  beg0: minimum value of slip and rake at each subfault
!  dp0: delta slip and delta rake at each subfault
!  np0: allowed amount of options for slip and rake values at each subfault
!  
   implicit none
   real :: time_min0(max_subfaults), time_max0(max_subfaults)
   real :: time_ref0(max_subfaults)
   real :: beg0(max_subfaults2), dp0(max_subfaults2)
   integer :: np0(max_subfaults2)
   time_min0(:) = time_min(:)
   time_max0(:) = time_max(:)
   time_ref0(:) = time_ref(:)
   beg0(:) = beg(:)
   dp0(:) = dp(:)
   np0(:) = np(:)
   end subroutine get_space

   
   subroutine get_borders(rake_min0, nleft0, nright0, nup0, ndown0)
!
!  Args:
!  rake_min0: minimum rake value allowed
!  nleft0: information about subfault to the left of current subfault
!  nright0: information about subfault to the right of current subfault
!  nup0: information about subfault above current subfault
!  ndown0: information about subfault below current subfault
!  
   implicit none
   integer :: nleft0(3, max_subfaults), nright0(3, max_subfaults), & 
   & nup0(3, max_subfaults), ndown0(3, max_subfaults)
   real :: rake_min0
   nleft0(:, :) = nleft(:, :)
   nright0(:, :) = nright(:, :)
   nup0(:, :) = nup(:, :)
   ndown0(:, :) = ndown(:, :)
   rake_min0 = rake_min
   end subroutine get_borders


   subroutine deallocate_ps()
   implicit none
   deallocate(point_sources)
   end subroutine deallocate_ps


end module model_parameters
