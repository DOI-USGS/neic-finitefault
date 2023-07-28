module save_forward


   use constants, only : max_seg, max_subf, wave_pts2, wave_pts, max_stations, n_data, twopi, &
           &      pi, max_rise_time_range, max_subfaults2, dpi 
   use wavelet_param, only : get_data_param 
   use get_stations_data, only : get_properties, llove, io_up
   use retrieve_gf, only : green_stk, green_dip
   use rise_time, only : source, fourier_asym_cosine, realtr, fft
   implicit none
   integer, parameter :: nnsta_tele = 80
   integer :: nxs_sub(max_seg), nys_sub(max_seg), nx_p, ny_p, msou, segments, subfaults
   real :: dta, ta0, dt_channel(max_stations)
   integer :: lnpt, max_freq, nlen
   character(len=15) :: sta_name(max_stations)
   character(len=3) :: component(max_stations)


contains


   subroutine saveforward_set_fault_parameters()
   use model_parameters, only : get_rise_time, get_segments, get_subfaults
   implicit none
   real :: dip(max_seg), strike(max_seg), delay_seg(max_seg), dxs, dys
   integer :: cum_subfaults(max_seg)
   real :: v_min, v_max, v_ref
   call get_rise_time(ta0, dta, msou)
   call get_segments(nxs_sub, nys_sub, dip, strike, delay_seg, segments, subfaults, cum_subfaults)
   call get_subfaults(dxs, dys, nx_p, ny_p, v_min, v_max, v_ref)
   end subroutine saveforward_set_fault_parameters


   subroutine saveforward_set_data_properties()
   implicit none
   integer :: jmin, jmax, channels
   call get_data_param(lnpt, jmin, jmax, nlen, max_freq)
   call get_properties(sta_name, component, dt_channel, channels) 
   end subroutine saveforward_set_data_properties


   subroutine write_forward(slip, rake, rupt_time, trise, tfall, &
      &  strong, cgps, body, surf, dart)
!
!  Args:
!  strong: True if strong motion data are used, False otherwise
!  cgps: True if cGPS data are used, False otherwise
!  body: True if body wave data are used, False otherwise
!  surf: True if surface wave data are used, False otherwise
!  dart: True if DART data are used, False otherwise
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  trise: array with model risetime values for all subfaults
!  tfall: array with model falltime values for all subfaults
!
!
!  Here, we write the forward solution given a kinematic model, for all specified 
!       data types.
!  
   implicit none
   integer first, last
   logical :: strong, cgps, body, surf, dart
   real slip(:), rake(:), rupt_time(:), tfall(:), trise(:)
   complex z0
!
   write(*,*)'Return synthetics from input kinematic model...'
   z0 = cmplx(0.0, 0.0)
   first = 0
   last = 0
!   write(*,*) dxs, dys
   if (strong) then
      call write_near_field_forward(slip, rake, rupt_time, trise, tfall, first, last, &
                                  & strong, .False.)
      first = last
   end if
   if (cgps) then
      call write_near_field_forward(slip, rake, rupt_time, trise, tfall, first, last, &
                                  & .False., cgps)
      first = last
   end if
   if (body) then
      call write_body_waves_forward(slip, rake, rupt_time, trise, tfall, first, last)
      first = last
   end if
   if (surf) then
      call write_surface_waves_forward(slip, rake, rupt_time, trise, tfall, first, last)
      first = last
   end if
   if (dart) then
      call write_dart_forward(slip, rake, rupt_time, trise, tfall, first, last)
      first = last
   end if
   end subroutine write_forward
   
   
   subroutine write_near_field_forward(slip, rake, rupt_time, trise, tfall, first, last, &
   &  strong, cgps)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  trise: array with model risetime values for all subfaults
!  tfall: array with model falltime values for all subfaults
!  first: number of initial channel
!  last: number of final channel
!  strong: True if strong motion data are used, False otherwise
!  cgps: True if cGPS data are used, False otherwise
!
   implicit none
   integer first, last, ll_g, isl, isr, channel_max, subfault, &
   &  jf, i, k, segment_subfault, segment, channel, n_chan, ixs, iys
   real slip(:), rake(:), rupt_time(:), &
   &  tfall(:), trise(:), real1(wave_pts2), imag1(wave_pts2), r, time, dt
   real*8 t1, t2, df
   complex forward(wave_pts2), z0, z
   complex :: source2(wave_pts, max_rise_time_range, max_rise_time_range)
   character(len=70) filename, filename2
   character(len=3) comp!component(max_stations), comp
   logical :: strong, cgps

   if (strong) write(*,*)'Return strong motion synthetics from input kinematic model...'
   if (cgps) write(*,*)'Return cGPS synthetics from input kinematic model...'
   z0 = cmplx(0.0, 0.0)
   filename = 'channels_strong.txt'
   if (cgps) filename = 'channels_cgps.txt'
   filename = trim(filename)
   open(9,file=filename,status='old')
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*) channel_max, n_chan
!   read(9,*)
!   do channel = 1, channel_max
!      read(9,*) int1, sta_name(channel), float1, float2, int2, component(channel), float3, int3
!   end do
   close(9)
!
!       make the rise time function
!   
   dt = dt_channel(first + 1)
   jf = 2**(lnpt-1)+1
   df = 1.d0/(2**lnpt)/dt
  
   do isl = 1, msou
      do isr = 1, msou
         source2(:, isl, isr) = cmplx(0.0, 0.0)
         t1 = ta0+(isl-1)*dta
         t2 = ta0+(isr-1)*dta
         if (t1 .lt. dt) t1 = dt
         if (t2 .lt. dt) t2 = dt
         do i = 1, jf
            call fourier_asym_cosine((i-1)*df, t1, t2, source2(i, isl, isr))
!                SOURCE(i, ISL, isr) = CMPLX(CR(i), CZ(i))*dt
         end do
      end do
   end do

!
!  end of rise time 
!       
   filename2 = 'synthetics_strong.txt'
   if (cgps) filename2 = 'synthetics_cgps.txt'
   filename2 = trim(filename2)
   open(18,file=filename2)
   k = 0
!       
!  set up the green function for every subfault
!  and calculate the initial value of objective function
!
   do channel = 1, n_chan
   
      ll_g = channel+first
      call create_waveform(slip, rake, rupt_time, trise, tfall, forward, source2, ll_g)
      comp = component(ll_g)

      do i = 1, jf
         if (i .le. max_freq) then
            real1(i) = real(forward(i))
            imag1(i) = aimag(forward(i))
         else
            real1(i) = 0.0
            imag1(i) = 0.0
         end if
      end do
       
      call realtr(real1, imag1, lnpt)
      call fft(real1, imag1, lnpt, 1.)
    
      write(18,*)nlen,dt,sta_name(ll_g),comp
      do k = 1, nlen
         write(18,*) real1(k), imag1(k)
      end do
   
   end do
   close(18)
   last = first+n_chan
   end subroutine write_near_field_forward


   subroutine write_body_waves_forward(slip, rake, rupt_time, trise, tfall, first, last)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  trise: array with model risetime values for all subfaults
!  tfall: array with model falltime values for all subfaults
!  first: number of initial channel
!  last: number of final channel
!
   implicit none
   integer n_chan, channel, ll_g, k, subfault, &
   &  segment, iys, jf, i, npxy, segment_subfault, ixs, isl, isr, nl, &
   &  first, last
   real slip(:), rake(:), rupt_time(:), tfall(:), trise(:), &
   &  dt, azim, w, real1(wave_pts2), imag1(wave_pts2), sinal, cosal
   real*8 t1, t2, df
   complex ::  z, z0, forward(wave_pts)
   complex :: source2(wave_pts, max_rise_time_range, max_rise_time_range)

   write(*,*)'Return body wave synthetics from input kinematic model...'
   open(9,file='channels_body.txt',status='old')
   read(9,*)
   read(9,*) 
   read(9,*) 
   read(9,*) n_chan
   close(9)

   z0 = cmplx(0.0, 0.0)

   dt = dt_channel(first + 1)
   jf = 2**(lnpt-1)+1
   df = 1.d0/(2**lnpt)/dt
   do isl = 1, msou
      do isr = 1, msou
         t1 = ta0+(isl-1)*dta
         t2 = ta0+(isr-1)*dta
         if (t1 .lt. dt) t1 = dt
         if (t2 .lt. dt) t2 = dt
         do i = 1, 2*max_freq
            call fourier_asym_cosine((i-1)*df, t1, t2, source2(i, isl, isr))
         end do
      end do
   end do
!
!  End of Rise Time 
!
   OPEN(18,FILE = 'synthetics_body.txt')

   k = 0
!       
!  Now, we compute the synthetic seismographs
!
   npxy = nx_p*ny_p
   do channel = 1, n_chan
      ll_g = first+channel
      call create_waveform(slip, rake, rupt_time, trise, tfall, forward, source2, ll_g)
      do i = 1, jf
         if (i .le. max_freq) then
            real1(i) = real(forward(i))
            imag1(i) = aimag(forward(i))
         else
            real1(i) = 0.0
            imag1(i) = 0.0
         end if
      end do
      call realtr(real1, imag1, lnpt)
      call fft(real1, imag1, lnpt, 1.0)
      nl = 2**lnpt
      if (llove(channel) .eq. 0) then
         write(18,*)nl,dt,sta_name(ll_g),'P'
      else
         write(18,*)nl,dt,sta_name(ll_g),'SH'
      end if
      do i = 1, nl
         write(18,*) real1(i), imag1(i)
      end do
   end do
!
   close(18)
   last = first+n_chan!nstaon
   end subroutine write_body_waves_forward


   subroutine write_surface_waves_forward(slip, rake, rupt_time, trise, tfall, first, last)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  trise: array with model risetime values for all subfaults
!  tfall: array with model falltime values for all subfaults
!  first: number of initial channel
!  last: number of final channel
!
   implicit none
   integer first, last, channel_max, subfault, &
   &  ll_g, isl, isr, jf, i, k, segment_subfault, segment, channel, n_chan, &
   &  iys, ixs!, io_up(max_stations)

   real slip(:), rake(:), rupt_time(:), tfall(:), trise(:), &
   &  real1(wave_pts2), imag1(wave_pts2), dt
   real*8 t1, t2, df

   complex z0, forward(wave_pts2), z
   complex :: source2(wave_pts, max_rise_time_range, max_rise_time_range)

   write(*,*)'Return long period surface wave synthetics from input kinematic model...'
   z0 = cmplx(0.0, 0.0)

   open(9,file='channels_surf.txt',status='old')
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*) channel_max, n_chan
   close(9)
!
! suppose the ni = u3e+11, then then moment of 1cm*1km^2 
! subfault is 3e+21. The gf is for Mo = 1e+20 dyne.cm
! the unit of surface wave is mm.
!
!       make the rise time function
!     
   dt = dt_channel(first + 1)
   jf = 2**(lnpt-1)+1
   df = 1.d0/(2**lnpt)/4.d0
   do isl = 1, msou
      do isr = 1, msou
         t1 = ta0+(isl-1)*dta
         t2 = ta0+(isr-1)*dta
         if (t1 .lt. dt) t1 = dt
         if (t2 .lt. dt) t2 = dt
         do i = 1, jf
            call fourier_asym_cosine((i-1)*df, t1, t2, source2(i, isl, isr))
         end do
      end do
   end do
!
!  end of rise time 
!       
   open(18,file='synthetics_surf.txt')
   k = 0
!       
!  set up the green function for every subfault
!  and calculate the initial value of objective function
!
   do channel = 1, n_chan
      ll_g = channel+first
      call create_waveform(slip, rake, rupt_time, trise, tfall, forward, source2, ll_g)
      do i = 1, jf
         if (i .le. max_freq) then
            real1(i) = real(forward(i))
            imag1(i) = aimag(forward(i))
         else
            real1(i) = 0.0
            imag1(i) = 0.0
         end if
      end do
     
      call realtr(real1, imag1, lnpt)
      call fft(real1, imag1, lnpt, 1.0)
   
      if (io_up(channel) .eq. 1) then
         write(18,*)nlen,dt,sta_name(ll_g),'P'
      else
         write(18,*)nlen,dt,sta_name(ll_g),'SH'
      end if
      do k = 1, nlen
         write(18,*) real1(k), imag1(k)
      end do
   end do   
   close(18)
   last = first+n_chan

   end subroutine write_surface_waves_forward


   subroutine write_dart_forward(slip, rake, rupt_time, trise, tfall, first, last)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  trise: array with model risetime values for all subfaults
!  tfall: array with model falltime values for all subfaults
!  first: number of initial channel
!  last: number of final channel
!
   implicit none
   integer first, last, ll_g, isl, isr, subfault, &
   &  jf, i, k, segment_subfault, segment, channel, n_chan, ixs, iys, channel_max
   real slip(:), rake(:), rupt_time(:), &
   &  tfall(:), trise(:), real1(wave_pts2), imag1(wave_pts2), &
   &  dt
   real*8 t1, t2, df
   complex forward(wave_pts2), z0, z
   complex :: source2(wave_pts, max_rise_time_range, max_rise_time_range)

   z0 = cmplx(0.0, 0.0)
   
   open(9,file='channels_dart.txt',status='old')
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*) channel_max, n_chan
   read(9,*)
!   do channel = 1, channel_max
!      read(9,*) int1, sta_name(channel), float1, float2, int2, component(channel), float3, int3
!   end do
   close(9)
!
!       make the rise time function
!     
   dt = dt_channel(first + 1)
   jf = 2**(lnpt-1)+1
   df = 1.d0/(2**lnpt)/dt
   do isl = 1, msou
      do isr = 1, msou
         t1 = ta0+(isl-1)*dta
         t2 = ta0+(isr-1)*dta
         if (t1 .lt. dt) t1 = dt
         if (t2 .lt. dt) t2 = dt
         do i = 1, jf
            call fourier_asym_cosine((i-1)*df, t1, t2, source2(i, isl, isr))
         end do
      end do
   end do
!
!  end of rise time 
!       
   open(18,file='synthetics_dart.txt')
   k = 0
!       
!  set up the green function for every subfault
!  and calculate the initial value of objective function
!
   do channel = 1, n_chan
      ll_g = channel+first
      call create_waveform(slip, rake, rupt_time, trise, tfall, forward, source2, ll_g)
   
      do i = 1, jf
         real1(i) = real(forward(i))
         imag1(i) = aimag(forward(i))
      end do
 
      call realtr(real1, imag1, lnpt)
      call fft(real1, imag1, lnpt, 1.)
   
      write(18,*)nlen,dt,sta_name(ll_g),'dart'
      do k = 1, nlen
         write(18,*) real1(k), imag1(k)
      end do
   end do   
   close(18)
   last = first+n_chan
   end subroutine write_dart_forward


   subroutine create_waveform(slip, rake, rupt_time, trise, tfall, forward, source2, ll_g)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  trise: array with model risetime values for all subfaults
!  tfall: array with model falltime values for all subfaults
!  forward: synthetic for current channel given the input model
!  source2: source function for given data type
!  ll_g: index of current channel
!
   implicit none
   integer ll_g, isl, isr, jf, segment, &
   &  i, k, subfault, ixs, iys, segment_subfault
   real slip(:), rake(:), rupt_time(:), &
   &  tfall(:), trise(:), a, b, ww, dt, rake2
   real*8 df
   complex :: forward(:), source2(:, :, :)
   complex :: z0, z
 
   z0 = cmplx(0.0, 0.0)
   dt = dt_channel(ll_g)
   jf = 2**(lnpt-1)+1
   df = 1.d0/(2**lnpt)/dt
   
   do i = 1, wave_pts
      forward(i) = z0
   end do
   do subfault = 1, subfaults
      isl = int((trise(subfault)-ta0)/dta+0.5)+1
      isr = int((tfall(subfault)-ta0)/dta+0.5)+1
      rake2 = rake(subfault)*dpi
      a = sin(rake2)*slip(subfault)
      b = cos(rake2)*slip(subfault)
      do i = 1, max_freq
         ww = -(i-1)*twopi*df*rupt_time(subfault)
         z = cmplx(cos(ww), sin(ww))
         forward(i) = forward(i) &
      & +(a*green_dip(i, ll_g, subfault)+b*green_stk(i, ll_g, subfault)) &
      & *source2(i, isl, isr)*z
      enddo
   end do
   end subroutine create_waveform


end module save_forward
        
