module model_parameters


   use constants, only : max_seg, max_stk_psources, max_dip_psources, max_stk_subfaults, &
               &     max_dip_subfaults, max_subf, max_subfaults, max_subfaults2, max_psources
   use modelling_inputs, only : get_weights_moment_end, get_many_events
   implicit none
   real :: slip0(max_subfaults), rake0(max_subfaults), rupt_time0(max_subfaults)
   real :: t_rise0(max_subfaults), t_fall0(max_subfaults)
   integer :: segments, nxs_sub(max_seg), nys_sub(max_seg), nx_p, ny_p, nxs0, nys0
   real :: dip(max_seg), strike(max_seg), delay_seg(max_seg), t_latest
   real, allocatable :: point_sources(:, :, :)
   real :: shear(max_subfaults)
   real :: dxs, dys, v_ref, v_min, v_max, lower_bound, upper_bound
   real :: rise_param1, rise_param2
   real :: time_min(max_subfaults), time_max(max_subfaults), rake_min
   real :: time_ref(max_subfaults)
   integer :: windows
   real :: minimum(max_subfaults2), delta(max_subfaults2)
   integer :: n_values(max_subfaults2)
   integer :: subfaults, cum_subfaults(max_seg), events
   logical :: segment_in_event(max_seg, 10)
   logical :: subfault_in_event(max_subfaults, 10)
!
! for regularization
!
   integer :: nleft(3, max_subfaults), nright(3, max_subfaults), & 
   & nup(3, max_subfaults), ndown(3, max_subfaults)


contains

   
   subroutine modelparameters_set_procedure_param()
   implicit none
   real :: real0, real1, real2, real3, real4
   call get_weights_moment_end(real0, real1, real2, real3, t_latest, real4)
   end subroutine modelparameters_set_procedure_param

   
   subroutine modelparameters_set_events()
   implicit none
   real :: moment_event(10)
   call get_many_events(events, moment_event)
   end subroutine modelparameters_set_events


   subroutine get_faults_data()
   implicit none
   integer int1, int2, int3, psource, subfaults_seg, &
   &  segment, i, j, nx_c, ny_c, subfault
   real dist, t_ref, t_max, t_min, delta, dip1, strike1, float1
   real :: shear2(max_subf, max_seg)
   allocate(point_sources(7, max_psources, max_subfaults))
!
!     Input Fault position to memory
!
   write(*,*)'Read and store fault segments data to memory...'
   open(12, file='fault&rise_time.txt', status='old')
   read(12,*) nxs0, nys0, float1
   read(12,*) segments, dxs, dys, nx_p, ny_p, v_min, v_max, lower_bound, upper_bound
   read(12,*) rise_param1, rise_param2, windows, v_ref, int1
   subfault = 0
   cum_subfaults(:) = 0
   do segment = 1, segments
      cum_subfaults(segment) = subfault
      read(12,*) int1, dip(segment), strike(segment)
      read(12,*) nxs_sub(segment), nys_sub(segment), delay_seg(segment)
      do i = 1, nxs_sub(segment) * nys_sub(segment)
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
      read(12,*) int1, dip1, strike1
      if ((abs(dip1-dip(segment)) .gt. 1.e-2).or. &
      &  (abs(strike1-strike(segment)) .gt. 1.e-2)) then
         write(*,*)'Fault mechanism in Fault.pos is not matched with that in Fault.das'
         write(*,*) segment
         write(*,*) dip1, dip(segment)
         write(*,*) strike1, strike(segment)
      end if

      subfaults_seg = nxs_sub(segment)*nys_sub(segment)
      do i = 1, subfaults_seg
         subfault = subfault+1
         do psource = 1, ny_p*nx_p
            read(12,*)(point_sources(j, psource, subfault), j = 1, 7)
         end do
      end do
   end do
   close(12)
!
!  minimum and maximum subfault rupture arrival time
!
   nx_c = max(int(nx_p/2.0+0.51), 1)
   ny_c = max(int(ny_p/2.0+0.51), 1)
   psource = (ny_c-1)*nx_p+nx_c
   do subfault = 1, subfaults
      dist = point_sources(4, psource, subfault)
      t_ref = dist/v_ref
      t_max = dist/v_min
      t_min = dist/v_max
      if (t_ref .gt. t_latest) t_ref = t_latest
      if (t_min .gt. t_latest) t_min = t_latest
      if (t_max .gt. t_latest) t_max = t_latest
      delta = t_min-t_ref
      if (lower_bound .lt. delta) then
         time_min(subfault) = delta
      else
         time_min(subfault) = lower_bound
      end if
!      time_ref(kxy, segment) = time_min(kxy, segment) + t_ref
      delta = t_max - t_ref
      if (upper_bound .gt. delta) then
         time_max(subfault) = delta
      else
         time_max(subfault) = upper_bound
      end if
   end do
!  
!  Input shear modulous Model
!
   open(12, file='shear_model.txt', status='old')
   read(12,*) int1
   if (int1.ne.segments) then
      write(*,*)'Amount of fault segments is different between mu file and Fault file'
      stop
   end if
   subfault = 0
   do segment = 1, segments
      read(12,*) int1, int2, int3
      subfaults_seg = nxs_sub(segment)*nys_sub(segment)
      read(12,*)(shear2(i, segment), i = 1, subfaults_seg)
      do i = 1, subfaults_seg
         subfault = subfault + 1
         shear(subfault) = shear2(i, segment) 
      enddo
   end do 
   close(12)
   end subroutine get_faults_data


   subroutine events_segments()
   implicit none
   integer :: event_segment(max_seg)
   integer :: integer1, integer2, i, k, segment, subfault0
   segment_in_event(:, :) = .False.
   subfault_in_event(:, :) = .False.
   open(20, file='segments_events.txt', status='old')
   read(20,*)
   do i=1, segments
      read(20, *) integer1, event_segment(i) 
   enddo
   close(20)
   do i = 1, events
      do segment = 1, segments
         subfault0 = cum_subfaults(segment)
         if (event_segment(segment) .eq. i) then
            segment_in_event(segment, i) = .True.
            do k = 1, nxs_sub(segment) * nys_sub(segment)
               subfault_in_event(k+subfault0, i) = .True.
            enddo
         endif
      enddo
   enddo
   end subroutine events_segments


   subroutine write_model(slip, rake, rupt_time, trise, tfall, use_waveforms)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  trise: array with model risetime values for all subfaults
!  tfall: array with model falltime values for all subfaults
!  use_waveforms: True if some waveform data used in modelling, False otherwise
!
   real :: slip(:), rake(:), rupt_time(:)
   real :: trise(:), tfall(:)
   real :: latitude_ep, longitude_ep, lat0, lon0, depth0, t_ref, moment_sol
   integer :: segment, iys, ixs, iy, ix, subfault, psource
   logical :: use_waveforms
   latitude_ep = 0.0
   longitude_ep = 0.0
   do subfault = 1, subfaults
      do psource = 1, ny_p*nx_p
!
! we have found the epicenter
!
         if (abs(point_sources(4, psource, subfault)) .le. 1e-3) then 
             latitude_ep = point_sources(1, psource, subfault)
             longitude_ep = point_sources(2, psource, subfault)
         end if
      end do
   end do


   open(13, file='Solucion.txt')
   write(13,*)'#Total number of fault_segments=', segments
   subfault = 0
   do segment = 1, segments
      write(13,131)'#Fault_segment =', segment, ' nx(Along-strike)=', &
     &  nxs_sub(segment), ' Dx = ', dxs, 'km ny(downdip)= ', nys_sub(segment), &
     & ' Dy = ', dys, 'km'
131           format(a, i4, a, i3, a, f5.2, a, i3, a, f5.2, a)
      write(13,132)'#Boundary of Fault_segment ', segment, &
     & '. EQ in cell (', nxs0, ',', nys0, '). Lon:', longitude_ep, &
     & '   Lat:', latitude_ep
132           format(a, i4, a, i2, a, i2, a, f10.4, a, f10.4)
      write(13,*)'#Lon.  Lat.  Depth'
      subfault = cum_subfaults(segment) + 1
      lat0 = point_sources(1, 1, subfault)
      lon0 = point_sources(2, 1, subfault)
      depth0 = point_sources(3, 1, subfault)
      write(13,*) lon0, lat0, depth0
      subfault = cum_subfaults(segment)
      subfault = subfault + (nys_sub(segment)-1)*nxs_sub(segment) + 1
      psource = (ny_p-1)*nx_p + 1
      lat0 = point_sources(1, psource, subfault)
      lon0 = point_sources(2, psource, subfault)
      depth0 = point_sources(3, psource, subfault)
      write(13,*) lon0, lat0, depth0 
      subfault = cum_subfaults(segment)
      subfault = subfault + nxs_sub(segment)*nys_sub(segment)
      psource = nx_p*ny_p
      lat0 = point_sources(1, psource, subfault)
      lon0 = point_sources(2, psource, subfault)
      depth0 = point_sources(3, psource, subfault)
      write(13,*) lon0, lat0, depth0
      subfault = cum_subfaults(segment) + nxs_sub(segment)
      psource = nx_p
      lat0 = point_sources(1, psource, subfault)
      lon0 = point_sources(2, psource, subfault)
      depth0 = point_sources(3, psource, subfault)
      write(13,*) lon0, lat0, depth0
      subfault = cum_subfaults(segment) + 1
      lat0 = point_sources(1, 1, subfault)
      lon0 = point_sources(2, 1, subfault)
      depth0 = point_sources(3, 1, subfault)
      write(13,*) lon0, lat0, depth0
      write(13,*)'#Lat. Lon. depth slip rake strike dip t_rup t_ris t_fal mo'
      ix = int(nx_p / 2) + 1
      iy = int(ny_p / 2) + 1
      psource = (iy-1)*nx_p + ix 
      subfault = cum_subfaults(segment)
      if (use_waveforms) then
         do iys = 1,nys_sub(segment)*nxs_sub(segment)
            subfault = subfault + 1
            t_ref = point_sources(5, psource, subfault)
            t_ref = min(t_ref, t_latest)
            moment_sol = slip(subfault) * shear(subfault) * dxs * dys * (1e10)
            lat0 = point_sources(1, 1, subfault)
            lon0 = point_sources(2, 1, subfault)
            depth0 = point_sources(3, 1, subfault)
            write(13, 133) lat0, lon0, depth0, &
         &  slip(subfault), rake(subfault), strike(segment), dip(segment), &
         &  rupt_time(subfault) + t_ref + delay_seg(segment), trise(subfault), &
         &  tfall(subfault), moment_sol
133  format(f14.6, f14.6, f14.6, f14.6, f14.6, f14.6, f14.6, f14.6, f14.6, f14.6, e14.6)
         end do
      else
         do iys = 1,nys_sub(segment)*nxs_sub(segment)
            subfault = subfault + 1
            t_ref = point_sources(5, psource, subfault)
            t_ref = min(t_ref, t_latest)
            moment_sol = slip(subfault) * shear(subfault) * dxs * dys * (1e10)
            lat0 = point_sources(1, 1, subfault)
            lon0 = point_sources(2, 1, subfault)
            depth0 = point_sources(3, 1, subfault)
            write(13, 133) lat0, lon0, depth0, &
         &  slip(subfault), rake(subfault), strike(segment), dip(segment), &
         &  0.0, 0.0, 0.0, moment_sol
         end do
      endif
   end do
   close(13)
   end subroutine write_model


   subroutine get_model_space()
   implicit none
!
!  boundary conditions? or what?
!
   integer :: nx0, ny0
   parameter(nx0 = max_stk_subfaults*2, ny0 = max_dip_subfaults*2)
   integer :: surf_disp, nblock, parameters, k, segment
   integer :: right_boundary, left_boundary, updip_boundary, downdip_boundary 
   real :: delt_x, delt_y, zmed_max, zmed_min, zleft_max, zleft_min
   real :: zright_max, zright_min, zup_max, zup_min, zdown_max, zdown_min
   real :: angle_max, angle_min, vel_max, vel_min
   real :: step1, step2
   integer :: int0, steps0, steps1, steps2, steps3, steps4, steps5, steps6, int1, int2
   integer :: nx, ny, i, j, segment2
   real :: surface(1000, 4), slip_boundary(nx0, ny0, 3), upper_bounds(5)
   real :: model_boundary(max_subfaults2, 4)
   open(17, file='model_space.txt', status='old')
   read(17,*) surf_disp
   if (surf_disp .eq. 1) then
      open(16, file='surface.constrain', status='old')
      read(16,*) nblock, surf_disp
      do i = 1, nblock
         read(16,*)(surface(i, k), k = 1, 4)
!
!  surface(i,1) is which segment is affected
!  surface(i,1) is which updip subfault in strike direction is affected
!  surface(i,2) is minimum allowed slip at said subfault
!  surface(i,3) is maximum allowed slip at said subfault
!
      end do
      close(16)
   end if
   parameters = 0
   k = 0
   do segment = 1, segments
      read(17,*) int0
      read(17,*) delt_x, delt_y
      read(17,*) left_boundary, right_boundary, updip_boundary, downdip_boundary
      read(17,*) zmed_max, zmed_min, steps0
      read(17,*) zleft_max, zleft_min, steps1
      read(17,*) zright_max, zright_min, steps2
      read(17,*) zup_max, zup_min, steps3
      read(17,*) zdown_max, zdown_min, steps4
      read(17,*) angle_max, angle_min, steps5
      read(17,*) vel_max, vel_min, steps6
      read(17,*) int1, int2
      nx = nxs_sub(segment)
      ny = nys_sub(segment)
      rake_min = angle_min
!
!  check whether the contrains is frighting each other
!
      if (right_boundary .eq. 1 .or. left_boundary .eq. 1) then
         step1 = (zmed_max-zleft_max)/(int((nx+1)/2))
         step2 = (zmed_max-zright_max)/(int((nx+1)/2))
         if (delt_x .lt. step1 .or. delt_x .lt. step2) then
            write(*,*)'the constrain of left or right boundary is error'
            stop
         end if
      end if
      if (downdip_boundary .eq. 1 .or. updip_boundary .eq. 1) then
         step1 = (zmed_max-zup_max)/(int((ny+1)/2))
         step2 = (zmed_max-zdown_max)/(int((ny+1)/2))
         if (delt_y .lt. step1 .or. delt_y .lt. step2) then
            write(*,*)'the constrain of left or right boundary is error'
            stop
         end if
      end if

!  First:  Give the range of mediate part of space

      do i = 2, nxs_sub(segment)-1
         do j = 2, nys_sub(segment)-1
            slip_boundary(i, j, 1) = zmed_min
            slip_boundary(i, j, 2) = zmed_max
            slip_boundary(i, j, 3) = steps0
         end do
      end do
!
!     Second: Give the value range of left and right part
!
      do j = 1, nys_sub(segment)
         if (left_boundary .eq. 1) then
            slip_boundary(1, j, 1) = zleft_min
            slip_boundary(1, j, 2) = zleft_max
            slip_boundary(1, j, 3) = steps1
         else
            slip_boundary(1, j, 1) = zmed_min
            slip_boundary(1, j, 2) = zmed_max
            slip_boundary(1, j, 3) = steps0
         end if
         if (right_boundary .eq. 1) then
            slip_boundary(nx, j, 1) = zright_min
            slip_boundary(nx, j, 2) = zright_max
            slip_boundary(nx, j, 3) = steps2
         else
            slip_boundary(nx, j, 1) = zmed_min
            slip_boundary(nx, j, 2) = zmed_max
            slip_boundary(nx, j, 3) = steps0
         end if
      end do
!
!  Give the value range of up and down part
!
      do i = 2, nxs_sub(segment)-1
         if (updip_boundary .eq. 1) then
            slip_boundary(i, 1, 1) = zup_min
            slip_boundary(i, 1, 2) = zup_max
            slip_boundary(i, 1, 3) = steps3
         else
            slip_boundary(i, 1, 1) = zmed_min
            slip_boundary(i, 1, 2) = zmed_max
            slip_boundary(i, 1, 3) = steps0
         end if
         if (downdip_boundary .eq. 1) then
            slip_boundary(i, ny, 1) = zdown_min
            slip_boundary(i, ny, 2) = zdown_max
            slip_boundary(i, ny, 3) = steps4
         else
            slip_boundary(i, ny, 1) = zmed_min
            slip_boundary(i, ny, 2) = zmed_max
            slip_boundary(i, ny, 3) = steps0
         end if
      end do
!
!  upper left corner
!
      if ((updip_boundary .eq. 1) .and. (left_boundary .eq. 1)) then
         slip_boundary(1, 1, 1) = min(zup_min, zleft_min)
         slip_boundary(1, 1, 2) = min(zup_max, zleft_max)
         slip_boundary(1, 1, 3) = min(steps3, steps1)
      elseif (updip_boundary .eq. 1) then
         slip_boundary(1, 1, 1) = min(zup_min, zmed_min) 
         slip_boundary(1, 1, 2) = min(zup_max, zmed_max)
         slip_boundary(1, 1, 3) = min(steps3, steps0)
      elseif (left_boundary .eq. 1) then
         slip_boundary(1, 1, 1) = min(zleft_min, zmed_min) 
         slip_boundary(1, 1, 2) = min(zleft_max, zmed_max)
         slip_boundary(1, 1, 3) = min(steps1, steps0)
      else
         slip_boundary(1, 1, 1) = zmed_min
         slip_boundary(1, 1, 2) = zmed_max
         slip_boundary(1, 1, 3) = steps0
      end if
!
!  upper right corner
!
      if ((updip_boundary .eq. 1) .and. (right_boundary .eq. 1)) then
         slip_boundary(nx, 1, 1) = min(zup_min, zright_min)
         slip_boundary(nx, 1, 2) = min(zup_max, zright_max)
         slip_boundary(nx, 1, 3) = min(steps3, steps2)
      elseif (updip_boundary .eq. 1) then
         slip_boundary(nx, 1, 1) = min(zup_min, zmed_min) 
         slip_boundary(nx, 1, 2) = min(zup_max, zmed_max)
         slip_boundary(nx, 1, 3) = min(steps3, steps0)
      elseif (right_boundary .eq. 1) then
         slip_boundary(nx, 1, 1) = min(zright_min, zmed_min) 
         slip_boundary(nx, 1, 2) = min(zright_max, zmed_max)
         slip_boundary(nx, 1, 3) = min(steps2, steps0)
      else
         slip_boundary(nx, 1, 1) = zmed_min
         slip_boundary(nx, 1, 2) = zmed_max
         slip_boundary(nx, 1, 3) = steps0
      end if
!
!  lower left corner
!
      if ((downdip_boundary .eq. 1) .and. (left_boundary .eq. 1)) then
         slip_boundary(1, ny, 1) = min(zdown_min, zleft_min)
         slip_boundary(1, ny, 2) = min(zdown_max, zleft_max)
         slip_boundary(1, ny, 3) = min(steps4, steps1)
      elseif (downdip_boundary .eq. 1) then
         slip_boundary(1, ny, 1) = min(zdown_min, zmed_min) 
         slip_boundary(1, ny, 2) = min(zdown_max, zmed_max)
         slip_boundary(1, ny, 3) = min(steps4, steps0)
      elseif (left_boundary .eq. 1) then
         slip_boundary(1, ny, 1) = min(zleft_min, zmed_min) 
         slip_boundary(1, ny, 2) = min(zleft_max, zmed_max)
         slip_boundary(1, ny, 3) = min(steps1, steps0)
      else
         slip_boundary(1, ny, 1) = zmed_min
         slip_boundary(1, ny, 2) = zmed_max
         slip_boundary(1, ny, 3) = steps0
      end if
!
!  lower right corner
!
      if ((downdip_boundary .eq. 1) .and. (right_boundary .eq. 1)) then
         slip_boundary(nx, ny, 1) = min(zdown_min, zright_min)
         slip_boundary(nx, ny, 2) = min(zdown_max, zright_max)
         slip_boundary(nx, ny, 3) = min(steps4, steps2)
      elseif (downdip_boundary .eq. 1) then
         slip_boundary(nx, ny, 1) = min(zdown_min, zmed_min) 
         slip_boundary(nx, ny, 2) = min(zdown_max, zmed_max)
         slip_boundary(nx, ny, 3) = min(steps4, steps0)
      elseif (right_boundary .eq. 1) then
         slip_boundary(nx, ny, 1) = min(zright_min, zmed_min) 
         slip_boundary(nx, ny, 2) = min(zright_max, zmed_max)
         slip_boundary(nx, ny, 3) = min(steps2, steps0)
      else
         slip_boundary(nx, ny, 1) = zmed_min
         slip_boundary(nx, ny, 2) = zmed_max
         slip_boundary(nx, ny, 3) = steps0
      end if
!
!  Recheck the range of mediate part
!
      do i = 2, nxs_sub(segment)-1
         do j = 2, nys_sub(segment)-1
            upper_bounds(1) = (i-1)*delt_x+zleft_max
            upper_bounds(2) = (nxs_sub(segment)-i)*delt_x+zright_max
            upper_bounds(3) = (j-1)*delt_y+zup_max
            upper_bounds(4) = (nys_sub(segment)-j)*delt_y+zdown_max
            upper_bounds(5) = slip_boundary(i, j, 2)
            slip_boundary(i, j, 2) = minval(upper_bounds)
         end do
      end do
!  
!  check the surface constrain
!  
      if (surf_disp .eq. 1) then
         do i = 1, nblock
            segment2 = int(surface(i, 1)+0.1)
            if (segment2 .eq. segment) then
               ny = 1
               nx = int(surface(i, 2)+0.1)
               slip_boundary(nx, ny, 1) = surface(i, 3)
               slip_boundary(nx, ny, 2) = surface(i, 4)
            end if
         end do
      end if
!
!  Change xy_range to slip_boundary    
!
      parameters = parameters+4*nxs_sub(segment)*nys_sub(segment)
      do j = 1, nys_sub(segment)
         do i = 1, nxs_sub(segment)
            k = k+1
            model_boundary(k, 1) = slip_boundary(i, j, 1)
            model_boundary(k, 2) = slip_boundary(i, j, 2)
            model_boundary(k, 3) = slip_boundary(i, j, 3)
            k = k+1
            model_boundary(k, 1) = angle_min
            model_boundary(k, 2) = angle_max
            model_boundary(k, 3) = steps5
            k = k+1
            model_boundary(k, 1) = vel_min
            model_boundary(k, 2) = vel_max
            model_boundary(k, 3) = steps6
            k = k+1
            model_boundary(k, 1) = 1
            model_boundary(k, 2) = windows
            model_boundary(k, 3) = windows
         end do
      end do
   end do
   close(17)

   do k = 1, parameters
      minimum(k) = model_boundary(k, 1)
      n_values(k) = int(model_boundary(k, 3) + 0.1)
      if (n_values(k) .gt. 1) then
         delta(k) = (model_boundary(k, 2)-model_boundary(k, 1))/(n_values(k)-1)
      else
         delta(k) = 0
      end if
   end do
   end subroutine get_model_space

   
   subroutine get_special_boundaries()
   implicit none
   integer :: special_subfaults, i, ixs, iys, segment
   integer :: subfault0, subfault1
!   real :: dd(max_subf, max_seg), aa(max_subf, max_seg)
!
!  special boundary
!
   open(12, file='special_model_space.txt', status='old')
   read(12,*) special_subfaults
   do i = 1, special_subfaults
      read(12,*) segment, ixs, iys
!      ll = 0
!      do segment = 1, segment-1
!         ll = ll+nxs_sub(segment)*nys_sub(segment)
!      end do
      subfault0 = cum_subfaults(segment)
      subfault1 = ixs+(iys-1)*nxs_sub(segment)
      subfault0 = subfault0+subfault1!ixs+(iys-1)*nxs_sub(segment)
      n_values(4*(subfault0-1)+1) = 2
      delta(4*(subfault0-1)+1) = 10
      minimum(4*(subfault0-1)+1) = 1
!      dd(subfault1, segment) = minimum(4*(subfault0-1)+1)
      n_values(4*(subfault0-1)+2) = 2
!      aa(subfault1, segment) = minimum(4*(subfault0-1)+2)
      n_values(4*(subfault0-1)+3) = 2
      n_values(4*(subfault0-1)+4) = 2
   end do
   close(12)
   end subroutine get_special_boundaries

   
   subroutine subfault_positions()
   implicit none
   integer k, special_subfaults
   integer segment2, subfault, subfault2
   integer ix, iy, stk_subf2, dip_subf2, change_neighbour
   character(len=80) char1
   integer segment
!
!  We detect adjacent subfaults and their relative location
!  
   subfault = 0
   do segment = 1, segments
      do iy = 1, nys_sub(segment)
         do ix = 1, nxs_sub(segment)
            subfault = subfault + 1
            nup(1, subfault) = segment
            nup(2, subfault) = ix
            nup(3, subfault) = iy-1
            ndown(1, subfault) = segment
            ndown(2, subfault) = ix
            ndown(3, subfault) = iy+1
            nleft(1, subfault) = segment
            nleft(2, subfault) = ix-1
            nleft(3, subfault) = iy
            nright(1, subfault) = segment
            nright(2, subfault) = ix+1
            nright(3, subfault) = iy
         end do
      end do
   end do
   open(22, file='regularization_borders.txt', status='old')
   do segment = 1, segments
      subfault = cum_subfaults(segment)
!      up 
      read(22,*)
      read(22,*) segment2, dip_subf2
      subfault2 = subfault
      do ix = 1, nxs_sub(segment)
         subfault2 = subfault2 + 1
         nup(1, subfault2) = segment2
         nup(2, subfault2) = ix
         nup(3, subfault2) = dip_subf2
      end do
!      down
      read(22,*) segment2, dip_subf2
      subfault2 = subfault + (nys_sub(segment)-1)*nxs_sub(segment)
      do ix = 1, nxs_sub(segment)
         subfault2 = subfault2 + 1
         ndown(1, subfault2) = segment2
         ndown(2, subfault2) = ix
         ndown(3, subfault2) = dip_subf2
      end do
!      left
      read(22,*) segment2, stk_subf2
      subfault2 = subfault + 1
      do iy = 1, nys_sub(segment)
         nleft(1, subfault2) = segment2
         nleft(2, subfault2) = stk_subf2
         nleft(3, subfault2) = iy
         subfault2 = subfault2 + nxs_sub(segment)
      end do
!      right
      read(22,*) segment2, stk_subf2
      subfault2 = subfault
      do iy = 1, nys_sub(segment)
         subfault2 = subfault2 + nxs_sub(segment)
         nright(1, subfault2) = segment2
         nright(2, subfault2) = stk_subf2
         nright(3, subfault2) = iy
      end do
   end do
   close(22)   
   open(22, file='special_regularization_borders.txt', status='old')
   read(22,*) special_subfaults
   do k = 1, special_subfaults
      read(22,'(a)')char1
! 1  format(a)
!      write(*,'(a)')aaaa
      read(22,*) segment, ix, iy
      read(22,*) change_neighbour, segment2, stk_subf2, dip_subf2
      subfault = cum_subfaults(segment) + nxs_sub(segment)*iy + ix
      if (change_neighbour .eq. 1) then
         nup(1, subfault) = segment2
         nup(2, subfault) = stk_subf2
         nup(3, subfault) = dip_subf2
      end if
      read(22,*) change_neighbour, segment2, stk_subf2, dip_subf2
      if (change_neighbour .eq. 1) then
         ndown(1, subfault) = segment2
         ndown(2, subfault) = stk_subf2
         ndown(3, subfault) = dip_subf2
      end if
      read(22,*) change_neighbour, segment2, stk_subf2, dip_subf2
      if (change_neighbour .eq. 1) then
         nleft(1, subfault) = segment2
         nleft(2, subfault) = stk_subf2
         nleft(3, subfault) = dip_subf2
      end if
      read(22,*) change_neighbour, segment2, stk_subf2, dip_subf2
      if (change_neighbour .eq. 1) then
         nright(1, subfault) = segment2
         nright(2, subfault) = stk_subf2
         nright(3, subfault) = dip_subf2
      end if
   end do
   close(22)
   end subroutine subfault_positions
   
   
   subroutine get_events_segments(segment_in_event0, subfault_in_event0)
!
!  Args:
!  segment_in_event0: True if given segment belongs to given event, False otherwise
!  subfault_in_event0: True if given subfault belongs to given event, False otherwise
!  
   implicit none
   logical :: segment_in_event0(max_seg, 10)
   logical :: subfault_in_event0(max_subfaults, 10)
   segment_in_event0(:, :) = segment_in_event(:, :)
   subfault_in_event0(:, :) = subfault_in_event(:, :)
   end subroutine get_events_segments


   subroutine get_rise_time(rise_param10, rise_param20, windows0)
!
!  Args:
!  rise_param10: first parameter of rise time function
!  rise_param20: second parameter of rise time function
!  windows: time windows for rise time function
!  
   implicit none
   real :: rise_param10, rise_param20
   integer :: windows0 
   rise_param10 = rise_param1
   rise_param20 = rise_param2
   windows0 = windows
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
   &  segments0, subfaults0, cum_subfaults0)
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
!  
   implicit none
   integer :: nxs_sub0(:), nys_sub0(:), subfaults0, segments0, cum_subfaults0(:)
   real :: dip0(:), strike0(:), delay_seg0(:)
   nxs_sub0(:) = nxs_sub(:)
   nys_sub0(:) = nys_sub(:)
   dip0(:) = dip(:)
   strike0(:) = strike(:)
   delay_seg0(:) = delay_seg(:)
   segments0 = segments
   subfaults0 = subfaults
   cum_subfaults0(:) = cum_subfaults(:)
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
   

   subroutine get_space(time_min0, time_max0, time_ref0, minimum0, delta0, n_values0)
!
!  Args:
!  time_min0: minimum allowed time of rupture initiation at each subfault
!  time_max0: maximum allowed time of rupture initiation at each subfault
!  time_ref0: reference time of rupture initiation at each subfault
!  minimum0: minimum value of slip and rake at each subfault
!  delta0: delta slip and delta rake at each subfault
!  n_valuesp0: allowed amount of options for slip and rake values at each subfault
!  
   implicit none
   real :: time_min0(max_subfaults), time_max0(max_subfaults)
   real :: time_ref0(max_subfaults)
   real :: minimum0(max_subfaults2), delta0(max_subfaults2)
   integer :: n_values0(max_subfaults2)
   time_min0(:) = time_min(:)
   time_max0(:) = time_max(:)
   time_ref0(:) = time_ref(:)
   minimum0(:) = minimum(:)
   delta0(:) = delta(:)
   n_values0(:) = n_values(:)
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
