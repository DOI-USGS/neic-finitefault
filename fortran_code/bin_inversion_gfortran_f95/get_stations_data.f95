module get_stations_data


   use constants, only : max_stations, max_seg, npuse, wave_pts2, wave_pts, n_data
   use wavelet_param, only : set_params
   use wavelets, only : wavelet_syn, cfft
   implicit none
   integer, parameter :: nnsta_tele = 150
   real :: weight(max_stations), dt_channel(max_stations)
   integer :: misfit_type(12, max_stations)
   integer :: t_min(max_stations), t_max(max_stations), channels
   real :: wavelet_weight(12, max_stations), observed(n_data, max_stations)
   character(len=15) :: sta_name(max_stations)
   character(len=3) :: component(max_stations)
   integer :: mmm(max_stations), llove(max_stations), io_up(max_stations), idata(max_stations), &
         &     disp_or_vel(max_stations)
   logical :: dart_channels(max_stations)
   logical :: cgnss_channels(max_stations)
   integer :: lnpt, nlen, jmin, jmax, max_freq, start(max_stations)
   integer :: event_sta(max_stations)


contains


   subroutine get_options(weight0, misfit_type0, &
   &  t_min0, t_max0, wavelet_weight0)
!
!  Args:
!  weight0: weight of channels in modelling
!  misfit_type0: type of objective function to be used
!  t_min0: start of waveform for each channel
!  t_max0: end of waveform for each channel
!  wavelet_weight0: weight of all wavelet scales for all channels
!
   implicit none
   real :: weight0(max_stations)
   integer :: misfit_type0(12, max_stations)
   integer :: t_min0(max_stations), t_max0(max_stations)
   real :: wavelet_weight0(12, max_stations)
   weight0(:) = weight(:)
   misfit_type0(:, :) = misfit_type(:, :)
   t_min0(:) = t_min(:)
   t_max0(:) = t_max(:)
   wavelet_weight0(:, :) = wavelet_weight(:, :)
   end subroutine get_options


   subroutine get_properties(sta_name0, component0, dt_channel0, channels0)
!
!  Args:
!  sta_name0: station name for each channel
!  component0: channel name for each channel
!  dt_channel0: sampling interval for each channel
!  channels0: amount of channels used in modelling
!
   implicit none
   character(len=15) :: sta_name0(max_stations)
   character(len=3) :: component0(max_stations)
   real :: dt_channel0(max_stations)
   integer :: channels0
   sta_name0(:) = sta_name(:)
   component0(:) = component(:)
   dt_channel0(:) = dt_channel(:)
   channels0 = channels
   end subroutine get_properties


   subroutine get_properties2(mmm0, llove0, idata0, disp_or_vel0)
!
!  Args:
!  mmm0: related to type of seismic source
!  llove0: whether P or S wave
!  idata0: ?
!  disp_or_vel0: whether to use displacement or velocity
!
   implicit none
   integer :: mmm0(max_stations), llove0(max_stations), idata0(max_stations)
   integer :: disp_or_vel0(max_stations)
   mmm0(:) = mmm(:)
   llove0(:) = llove(:)
   idata0(:) = idata(:)
   disp_or_vel0(:) = disp_or_vel(:)
   end subroutine get_properties2
   

   subroutine get_event_sta(event_sta0)
!
!  Args:
!  event_sta0: earthquake at which each channel belongs
!
   implicit none
   integer :: event_sta0(max_stations)
   event_sta0(:) = event_sta(:)
   end subroutine get_event_sta


   subroutine get_data(strong, cgnss, body, surf, dart)
!
!  Args:
!  strong: True if strong motion data are used, False otherwise
!  cgnss: True if cGNSS data are used, False otherwise
!  body: True if body wave data are used, False otherwise
!  surf: True if surface wave data are used, False otherwise
!  dart: True if DART data are used, False otherwise
!
   implicit none
!
!  Here, we load into memory, wavelet transform of observed data, and 
!       other properties of stations
!
   integer :: first, last
   logical :: strong, cgnss, body, surf, dart

   write(*,*)'Get stations metadata and waveforms and store them in memory...'
   first = 0
   last = 0
   dt_channel(:) = 0.0
   if (strong) then
      call get_near_field_stations(first, last, strong, .False.)
      first = last
   end if
   if (cgnss) then
      call get_near_field_stations(first, last, .False., cgnss)
      first = last
   end if
   if (body) then
      call get_body_waves_stations(first, last)
      first = last
   end if
   if (surf) then
      call get_surface_waves_stations(first, last)
      first = last
   end if
   if (dart) then
      call get_dart_stations(first, last)
      first = last
   end if
   channels = last
   end subroutine get_data


   subroutine get_near_field_stations(first, last, strong, cgnss)
!
!  Args:
!  first: number of initial channel
!  last: number of final channel
!  strong: True if strong motion data are used, False otherwise
!  cgnss: True if cGNSS data are used, False otherwise
!
   implicit none
   integer first, last, channel, channel0, misfit_type0(11), k, i, int0, &
   &  n_wave_weight, stations, channels0, int1, int2, used_channels
   real dt, weig(max_stations), wavelet_weight0(11), dt_sample, lat_s, lon_s
   logical, parameter :: dart = .False.
   character(len=20) filename, string2, string1
   character(len=30) event_file, channels_file, wavelets_file, waveforms_file
   logical :: is_file, strong, cgnss
!   character(len=6) sta_name(max_stations)
!   character(len=3) component(max_stations)

   if (strong) write(*,*)'Get strong motion stations metadata and waveforms...'
   if (cgnss) write(*,*)'Get cGNSS stations metadata and waveforms...'
   used_channels = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!

   channels_file = 'channels_strong.txt'
   if (cgnss) channels_file = 'channels_cgnss.txt'
   channels_file = trim(channels_file)
   open(9,file=channels_file,status='old')
   wavelets_file = 'wavelets_strong.txt'
   if (cgnss) wavelets_file = 'wavelets_cgnss.txt'
   wavelets_file = trim(wavelets_file)
   open(15,file=wavelets_file, status='old')
   read(15,*) jmin, jmax, max_freq
   read(15,*)
   read(15,*) n_wave_weight

   read(9,*)
   read(9,*) 
   read(9,*) lnpt, dt_sample
   read(9,*)
   nlen = 2**lnpt
   call set_params(lnpt, jmin, jmax, nlen, max_freq)
   dt = dt_sample

   read(9,*) stations, channels0
   read(9,*)
   call error1(first, channels0)
   call error2(n_wave_weight, channels0)
   channel0 = 0
 
   do i = 1, stations
      channel = i+first
      cgnss_channels(channel) = .False.
      if (cgnss) cgnss_channels(channel) = .True.
      dart_channels(channel) = .False.
      read(9,*) int0, sta_name(channel), lat_s, lon_s, int1, component(channel), weig(i), int2
      if (weig(i) .gt. 0) used_channels = used_channels + 1
   end do
   close(9)
!   write(*,*)'channels0: ', used_channels

   filename = 'waveforms_strong.txt'
   if (cgnss) filename = 'waveforms_cgnss.txt'
   filename = trim(filename)
   call get_waveforms(filename, channels0, first)

   do i = 1, stations

      channel0 = channel0+1
      channel = channel0+first
      weight(channel) = weig(i) / used_channels
      dt_channel(channel) = dt_sample
      
      do k = 1, 11
         misfit_type(k, channel) = 0
         wavelet_weight(k, channel) = 0
      end do
      read(15,*)
      read(15,*)(misfit_type0(k), k = 1, jmax)
      read(15,*)(wavelet_weight0(k), k = 1, jmax)
      do k = jmin, jmax
         misfit_type(k, channel) = misfit_type0(k)
         wavelet_weight(k, channel) = wavelet_weight0(k)
      end do
   end do
   close(15) 

   event_file = 'strong_motion_events.txt'
   if (cgnss) event_file = 'cgnss_events.txt'
   inquire(file = event_file, exist = is_file)
   if (is_file) then
      open(12, file=event_file, status='old')
      do i=1, channels0
         channel = i+first
         read(12,*)string1, string2, event_sta(channel)
      enddo
      close(12)    
   endif
 
   last = first + channels0
   end subroutine get_near_field_stations
   
   
   subroutine get_body_waves_stations(first, last)
!
!  Args:
!  first: number of initial channel
!  last: number of final channel
!
   implicit none
   integer first, last, channel, misfit_type0(11), k, i, n_wave_weight, channels0, &
   &  int0, used_channels, int1, love, int2
   real lat_sta, lon_sta, wavelet_weight0(11), dt, & 
   &  rang, az, angle, float1, float2
   logical, parameter :: cgnss=.False., dart = .False.
   character(len=20) filename, string1, string2
   character(len=30) event_file
   character(len=6)  earth, sttyp 
   character(len=14) fname
   logical :: is_file

   write(*,*)'Get body waves stations metadata and waveforms...'
   used_channels = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   open(9, file='channels_body.txt', status='old')
   open(14, file='body_wave_weight.txt', status='old')
   if (max_freq .gt. npuse) then
      write(*,*)'You should stop and check dimension sww,sws'
      stop
   end if
   read(9,*)
   read(9,*) int0, lnpt, dt
   read(9,*)
   read(9,*) channels0
   call error1(first, channels0)
   nlen = 2 ** lnpt

   open(15, file='wavelets_body.txt', status='old')
   read(15,*) jmin, jmax, max_freq
   read(15,*)
   read(15,*) n_wave_weight
   do i = 1, channels0
      channel = first+i
      read(14,*) weight(channel)
      if (weight(channel) .gt. 0) used_channels = used_channels + 1
   end do
   rewind 14
   call set_params(lnpt, jmin, jmax, nlen, max_freq)
!   write(*,*)'channels0: ', used_channels
   call error2(n_wave_weight, channels0)
   
   filename = 'waveforms_body.txt'
   filename = trim(filename)
   call get_waveforms(filename, channels0, first)!, cgnss, dart)

   do i = 1, channels0
      channel = i+first
      cgnss_channels(channel) = .False.
      dart_channels(channel) = .False.
      read(9,*) int1, earth, sttyp, sta_name(channel), fname, &
      & rang, az, lat_sta, lon_sta, angle, float1, &
      & mmm(i), disp_or_vel(i), float2, llove(i), int2, idata(i)
      if (idata(i) .gt. 0 .or. mmm(i) .eq. 3) cycle
      dt_channel(channel) = dt
      read(14,*) weight(channel)
      weight(channel) = weight(channel) / used_channels
      love = llove(i)
      component(channel) = 'P'
      if (love .gt. 0) component(channel) = 'SH'
      
      do k = 1, 11
         misfit_type(k, channel) = 0
         wavelet_weight(k, channel) = 0
      end do
      read(15,*)
      read(15,*)(misfit_type0(k), k = 1, jmax)
      read(15,*)(wavelet_weight0(k), k = 1, jmax)
      do k = jmin, jmax
         misfit_type(k, channel) = misfit_type0(k)
         wavelet_weight(k, channel) = wavelet_weight0(k)
      end do
   end do
   close(9)
   close(14)
   close(15)
   
   event_file = 'tele_events.txt'
   inquire(file = event_file, exist = is_file)
   if (is_file) then
      open(12, file=event_file, status='old')
      do i=1, channels0
         channel = i+first
         read(12,*)string1, string2, event_sta(channel)
      enddo
      close(12)    
   endif
   
   last = first + channels0
   end subroutine get_body_waves_stations


   subroutine get_surface_waves_stations(first, last)
!
!  Args:
!  first: number of initial channel
!  last: number of final channel
!
   implicit none
   integer first, last, channel, channel0, misfit_type0(11), k, i, int0, &
   &  n_wave_weight, stations, channels0, int1, int2, used_channels, &
   &  io_ns(max_stations), io_ew(max_stations), int3
   real lat_s, lon_s, dt, ang_ns, ang_ew, &
   &  weig(max_stations, 3), wavelet_weight0(11), dt_sample, &
   &  dip, rake, theta
   logical, parameter :: cgnss=.False., dart=.False.
   character(len=20) filename, string1, string2
   character(len=30) event_file
!   character(len=6) sta_name(max_stations)
   character(len=250) modes
   logical :: is_file

   write(*,*)'Get stations metadata and waveforms for long period surface waves...'
   used_channels = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
 
   open(9, file='channels_surf.txt', status='old')

   open(15, file='wavelets_surf.txt', status='old')
   read(15,*) jmin, jmax, max_freq
   read(15,'(a)')modes
   read(15,*) n_wave_weight
 
   read(9,*)
   read(9,*)
   read(9,*) dip, rake, theta, lnpt, dt_sample
   read(9,*)
   if (lnpt.ne.10) then
      write(*,*)"please check input LNPT"
   end if
   nlen = 2**lnpt
   dt = dt_sample
   call set_params(lnpt, jmin, jmax, nlen, max_freq)
  
   read(9,*) stations, channels0
   read(9,*)
   call error1(first, channels0)
   call error2(n_wave_weight, channels0)
!       
!       Here we read the green functions of long period surface waves
!  
   do i = 1, stations    
      channel = i+first
      cgnss_channels(channel) = .False.
      dart_channels(channel) = .False.
      read(9,*) int0, sta_name(channel), lat_s, lon_s, int1, &
      & io_up(i), io_ns(i), io_ew(i), ang_ns, ang_ew, &
      & int2,(weig(i, k), k = 1, 3), int3
      if ((weig(i, 1) .gt. 0) .and. (io_up(i) .eq. 1)) used_channels = used_channels + 1
      if ((weig(i, 2) .gt. 0) .and. (io_ns(i) .eq. 1)) used_channels = used_channels + 1
   end do
   close(9)
   
   filename = 'waveforms_surf.txt'
   filename = trim(filename)
   call get_waveforms(filename, channels0, first) 

   channel0 = 0
   do i = 1, stations
      channel0 = channel0+1
      channel = channel0+first
      dt_channel(channel) = dt_sample
      if (io_up(i) .eq. 1) then
         weight(channel) = weig(i, 1) / used_channels
         component(channel) = 'R'
      end if
      if (io_ns(i) .eq. 1) then
         weight(channel) = weig(i, 2) / used_channels
         component(channel) = 'L'
      end if
      
      do k = 1, 11
         misfit_type(k, channel) = 0
         wavelet_weight(k, channel) = 0
      end do
      read(15,*)
      read(15,*)(misfit_type0(k), k = 1, jmax)
      read(15,*)(wavelet_weight0(k), k = 1, jmax)
      do k = jmin, jmax
         misfit_type(k, channel) = misfit_type0(k)
         wavelet_weight(k, channel) = wavelet_weight0(k)
      end do
   end do
   close(15)   
   
   event_file = 'surf_events.txt'
   inquire(file = event_file, exist = is_file)
   if (is_file) then
      open(12, file=event_file, status='old')
      do i=1, channels0
         channel = i+first
         read(12,*)string1, string2, event_sta(channel)
      enddo
      close(12)    
   endif
   
   last = first+channels0
   end subroutine get_surface_waves_stations


   subroutine get_dart_stations(first, last)
!
!  Args:
!  first: number of initial channel
!  last: number of final channel
!
   implicit none
   integer first, last, channel, misfit_type0(11), k, i, int0, &
   &  n_wave_weight, stations, channels0, int1, int2, used_channels
   real :: lat_s, lon_s, dt, weig(max_stations), wavelet_weight0(11), dt_sample
   logical, parameter :: cgnss=.False., dart=.True.
   character(len=20) filename
!   character(len=6) sta_name(max_stations)
!   character(len=3) component(max_stations)

   used_channels = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   open(9, file='channels_dart.txt', status='old')
   open(15, file='wavelets_dart.txt', status='old')
   read(15,*) jmin, jmax, max_freq
   read(15,*)
   read(15,*) n_wave_weight

   read(9,*)
   read(9,*)
   read(9,*) lnpt, dt_sample
   read(9,*)
   nlen = 2**lnpt
   call set_params(lnpt, jmin, jmax, nlen, max_freq)
   dt = dt_sample

   read(9,*) stations, channels0
   read(9,*)
   call error1(first, channels0)
   call error2(n_wave_weight, channels0)
 
   do i = 1, stations
      channel = i+first
      cgnss_channels(channel) = .False.
      dart_channels(channel) = .True.
      read(9,*) int0, sta_name(channel), lat_s, lon_s, int1, component(channel), weig(i), int2
      if (weig(i) .gt. 0) used_channels = used_channels + 1
   end do
   close(9)
   
   filename = 'waveforms_dart.txt'
   filename = trim(filename)
   call get_waveforms(filename, stations, first) 

   do i = 1, stations

      channel = i+first
      weight(channel) = weig(i) / used_channels
      dt_channel(channel) = dt_sample
      
      do k = 1, 11
         misfit_type(k, channel) = 0
         wavelet_weight(k, channel) = 0
      end do
      read(15,*)
      read(15,*)(misfit_type0(k), k = 1, jmax)
      read(15,*)(wavelet_weight0(k), k = 1, jmax)
      do k = jmin, jmax
         misfit_type(k, channel) = misfit_type0(k)
         wavelet_weight(k, channel) = wavelet_weight0(k)
      end do
   end do
   close(15)
   last = first + channels0
   end subroutine get_dart_stations


   subroutine get_waveforms(filename, channels0, first)
!
!  Args:
!  filename: name of file with waveforms
!  channels0: amount of channels to be read
!  first: number of initial channel
!
   implicit none
   integer first, channel, i, j, channels0
   real :: dto
   character(len=40) string
   character(len=20) filename
  
   filename = trim(filename)
   open(13, file=filename, status='old')
   do i = 1, channels0
      channel = i+first
      t_min(channel) = 0
      read(13,*)
      read(13,*)
      read(13,*) string, dto
      read(13,*) string, start(channel)
      read(13,*) string, t_max(channel)
      read(13,*)
      read(13,*)(observed(j, channel), j = 1, start(channel))
   enddo
   close(13)
   end subroutine get_waveforms


   subroutine get_wavelet_obs(wave_obs, max_coeff, t_max_val)
!
!  Args:
!  wave_obs: wavelet coefficients for channel waveform
!  max_coeff: maximum wavelet coefficient
!  t_max_val: wavelet atom belonging to maximum coefficient
!
   implicit none
   real, intent(out) :: wave_obs(wave_pts2, max_stations), max_coeff(max_stations)
   integer, intent(out) :: t_max_val(max_stations)
   integer atom_max0, i, j, channel, start1, length, n_begin, n_delt, index0
   real real1(wave_pts2), imag1(wave_pts2), observed2(n_data), &
   &  max_coeff0, mean
   real :: coeffs_obs(n_data)
   logical :: cgnss, dart
!
! TODO : how to implement this in a more elegant way?
!
   do channel = 1, channels
      start1 = start(channel)
      cgnss = cgnss_channels(channel)
      dart = dart_channels(channel)
      mean = sum(observed(start1 - 20:start1, channel)) / 20.0
      index0 = min(start1, nlen - 80)
      do i = 1, wave_pts2
         real1(i) = 0.0 
         imag1(i) = 0.0
         if (cgnss) then
            if (i .lt. index0) then
               observed2(i) = observed(i, channel) 
            else
               if (i .le. nlen - 80) observed2(i) = mean
               if (i .gt. nlen - 80) observed2(i) = 0.0
            endif
         else
            if (i .lt. start1) then
               observed2(i) = observed(i, channel) 
            else
               observed2(i) = 0.0 
            endif
         endif
      end do
      do i = 1,nlen
         real1(i) = observed2(i)*nlen
         imag1(i) = 0.0
      enddo

      call cfft(real1, imag1, lnpt)
      call wavelet_syn(real1, imag1, coeffs_obs)

      observed2(max_freq:) = 0.0
      max_coeff0 = 0.0
      atom_max0 = 1
      do i = 1, nlen
         if (max_coeff0 .lt. abs(coeffs_obs(i))) then
            max_coeff0 = abs(coeffs_obs(i))
            atom_max0 = i
         end if
      end do
      if (dart .or. cgnss) then
         max_coeff0 = 0.0
         do j = jmin, jmax
            n_begin = 2**(j-1)
            n_delt = nlen/n_begin
            length = int(t_max(channel)/n_delt+0.5)-1
            do i =1, length
               if (max_coeff0 .lt. abs(coeffs_obs(n_begin+i))) then
                  max_coeff0 = abs(coeffs_obs(n_begin+i))
                  atom_max0 = n_begin+i
               endif
            enddo
         enddo
      endif

      max_coeff(channel) = max_coeff0
      t_max_val(channel) = atom_max0
      do i = 1, nlen
         wave_obs(i, channel) = coeffs_obs(i)/(max_coeff0)
      end do
   enddo
   end subroutine get_wavelet_obs 


   subroutine count_wavelets(used_data)
!
!  Args:
!  used_data: amount of wavelet coefficients used
!
   implicit none
   integer :: used_data
   integer :: j, channel, n_begin, n_delt, length
   do channel = 1, channels
      if (weight(channel) .gt. (1e-3)) then
         do j = jmin, jmax
            n_begin = 2**(j-1)
            n_delt = nlen/n_begin
            length = int(t_max(channel)/n_delt+0.5)-1
            if (wavelet_weight(j, channel) .gt. (1e-3)) then
               used_data = used_data + length
            endif
         end do
      endif
   end do
   end subroutine count_wavelets


   subroutine error1(first, channels0)
!
!  Args:
!  first: number of initial channel
!  channels0: amount of channels used by specific data type
!
   implicit none
   integer :: first, channels0
   if (channels0 + first .gt. max_stations) then
      write(*,*)'Error: Maximum allowed number of channels: ', max_stations
      write(*,*)'But amount of channels used for modelling is at least', channels0 + first
      write(*,*)'Please check maximum allowed number of channels'
      stop
   end if
   end subroutine error1
   

   subroutine error2(n_wave_weight, channels0)
!
!  Args:
!  n_wave_weight: amount of channels specified in wavelet file 
!  channels0: amount of channels used by specific data type
!
   implicit none
   integer :: n_wave_weight, channels0
   if (n_wave_weight.ne.channels0) then
      write(0,*)'n_wave_chan channels0',n_wave_weight,channels0
      write(0,*)'Mismatch in amount of stations betweeen wavelet file and stations file'
      stop
   end if
   end subroutine error2


end module get_stations_data
