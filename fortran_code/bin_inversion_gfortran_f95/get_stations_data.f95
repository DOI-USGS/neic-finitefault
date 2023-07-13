module get_stations_data


   use constants, only : max_stations, max_seg, npuse, wave_pts2, wave_pts, n_data
   use wavelet_param, only : set_params
   use wavelets, only : wavelet_obs
   implicit none
   integer, parameter :: nnsta_tele = 150
   real :: weight(max_stations), dt_channel(max_stations)
   integer :: misfit_type(12, max_stations)
   integer :: t_min(max_stations), t_max(max_stations), channels
   real :: wavelet_weight(12, max_stations), obse(n_data, max_stations)
   character(len=15) :: sta_name(max_stations)
   character(len=3) :: component(max_stations)
   integer :: mmm(max_stations), llove(max_stations), io_up(max_stations), idata(max_stations), &
         &     disp_or_vel(max_stations)
   logical :: dart_channels(max_stations)
   logical :: cgps_channels(max_stations)
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

   
   subroutine get_event_sta(event_sta0)
!
!  Args:
!  event_sta0: earthquake at which each channel belongs
!
   implicit none
   integer :: event_sta0(max_stations)
   event_sta0(:) = event_sta(:)
   end subroutine get_event_sta

   subroutine get_data(strong, cgps, body, surf, dart)
!
!  Args:
!  strong: True if strong motion data are used, False otherwise
!  cgps: True if cGPS data are used, False otherwise
!  body: True if body wave data are used, False otherwise
!  surf: True if surface wave data are used, False otherwise
!  dart: True if DART data are used, False otherwise
!
   implicit none
!
!  Here, we load into memory, wavelet transform of observed data, and 
!       other properties of stations
!
   integer :: ll_in, ll_out
   logical :: strong, cgps, body, surf, dart

   write(*,*)'Get stations metadata and waveforms and store them in memory...'
   ll_in = 0
   ll_out = 0
   dt_channel(:) = 0.0
   if (strong) then
      call get_near_field_stations(ll_in, ll_out, strong, .False.)
      ll_in = ll_out
   end if
   if (cgps) then
      call get_near_field_stations(ll_in, ll_out, .False., cgps)
      ll_in = ll_out
   end if
   if (body) then
      call get_body_waves_stations(ll_in, ll_out)
      ll_in = ll_out
   end if
   if (surf) then
      call get_surface_waves_stations(ll_in, ll_out)
      ll_in = ll_out
   end if
   if (dart) then
      call get_dart_stations(ll_in, ll_out)
      ll_in = ll_out
   end if
   channels = ll_out
   end subroutine get_data


   subroutine get_near_field_stations(ll_in, ll_out, strong, cgps)
!
!  Args:
!  ll_in: number of initial channel
!  ll_out: number of final channel
!  strong: True if strong motion data are used, False otherwise
!  cgps: True if cGPS data are used, False otherwise
!
   implicit none
   integer ll_in, ll_out, ll_g, io_chan, j_con(11), k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, &
   &  nos(max_stations), io_mod(max_stations), n_chan3
   real dt, weig(max_stations), j_wig(11), dt_sample, lat_s, lon_s
   logical, parameter :: dart = .False.
   character(len=20) filename, string2, string1
   character(len=30) event_file, channels_file, wavelets_file, waveforms_file
   logical :: is_file, strong, cgps
!   character(len=6) sta_name(max_stations)
!   character(len=3) component(max_stations)

   if (strong) write(*,*)'Get strong motion stations metadata and waveforms...'
   if (cgps) write(*,*)'Get cGPS stations metadata and waveforms...'
   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!

   channels_file = 'channels_strong.txt'
   if (cgps) channels_file = 'channels_cgps.txt'
   channels_file = trim(channels_file)
   open(9,file=channels_file,status='old')
   wavelets_file = 'wavelets_strong.txt'
   if (cgps) wavelets_file = 'wavelets_cgps.txt'
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

   read(9,*) ir_max, n_chan
   read(9,*)
   call error1(ll_in, n_chan)
   call error2(n_wave_weight, n_chan)
   ir = 0
   io_chan = 0
 
   do ir = 1, ir_max
      ll_g = ir+ll_in
      cgps_channels(ll_g) = .False.
      if (cgps) cgps_channels(ll_g) = .True.
      dart_channels(ll_g) = .False.
      read(9,*) no, sta_name(ll_g), lat_s, lon_s, io_mod(ir), component(ll_g), weig(ir), nos(ir)
      if (weig(ir) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   close(9)
!   write(*,*)'n_chan: ', n_chan3

   filename = 'waveforms_strong.txt'
   if (cgps) filename = 'waveforms_cgps.txt'
   filename = trim(filename)
   call get_waveforms(filename, ir_max, ll_in, ll_out)

   do ir = 1, ir_max

      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      weight(ll_g) = weig(ir) / n_chan3
      dt_channel(ll_g) = dt_sample
      
      do k = 1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
      read(15,*)
      read(15,*)(j_con(k), k = 1, jmax)
      read(15,*)(j_wig(k), k = 1, jmax)
      do k = jmin, jmax
         misfit_type(k, ll_g) = j_con(k)
         wavelet_weight(k, ll_g) = real(j_wig(k))
      end do
      do k = jmax+1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
   end do
   close(15) 

   event_file = 'strong_motion_events.txt'
   if (cgps) event_file = 'cgps_events.txt'
   inquire(file = event_file, exist = is_file)
   if (is_file) then
      open(12, file=event_file, status='old')
      do ir=1,ir_max
         ll_g = ir+ll_in
         read(12,*)string1, string2, event_sta(ll_g)
      enddo
      close(12)    
   endif
 
   ll_out = ll_in + n_chan
   end subroutine get_near_field_stations
   
   
   subroutine get_body_waves_stations(ll_in, ll_out)
!
!  Args:
!  ll_in: number of initial channel
!  ll_out: number of final channel
!
   implicit none
   integer ll_in, ll_out, ll_g, j_con(11), k, ir, n_wave_weight, n_chan, &
   &  nos, n_chan3, idts, love, int1 
   real lat_sta, lon_sta, j_wig(11), dt, & 
   &  rang, az, earth_angle, float1, float2
   logical, parameter :: cgps=.False., dart = .False.
   character(len=20) filename, string1, string2
   character(len=30) event_file
   character(len=6)  earth, sttyp 
   character(len=14) fname
   logical :: is_file

   write(*,*)'Get body waves stations metadata and waveforms...'
   n_chan3 = 0
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
   read(9,*) idts, lnpt, dt
   read(9,*)
   read(9,*) n_chan
   call error1(ll_in, n_chan)
   nlen = 2 ** lnpt

   open(15, file='wavelets_body.txt', status='old')
   read(15,*) jmin, jmax, max_freq
   read(15,*)
   read(15,*) n_wave_weight
   do ir = 1, n_chan
      ll_g = ll_in+ir
      read(14,*) weight(ll_g)
      if (weight(ll_g) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   rewind 14
   call set_params(lnpt, jmin, jmax, nlen, max_freq)
!   write(*,*)'n_chan: ', n_chan3
   call error2(n_wave_weight, n_chan)
   
   filename = 'waveforms_body.txt'
   filename = trim(filename)
   call get_waveforms(filename, n_chan, ll_in, ll_out)!, cgps, dart)

   do ir = 1, n_chan
      ll_g = ir+ll_in
      cgps_channels(ll_g) = .False.
      dart_channels(ll_g) = .False.
      read(9,*) nos, earth, sttyp, sta_name(ll_g), fname, &
      & rang, az, lat_sta, lon_sta, earth_angle, float1, &
      & mmm(ir), disp_or_vel(ir), float2, llove(ir), int1, idata(ir)
      if (idata(ir) .gt. 0 .or. mmm(ir) .eq. 3) cycle
      dt_channel(ll_g) = dt
      read(14,*) weight(ll_g)
      weight(ll_g) = weight(ll_g) / n_chan3
      love = llove(ir)
      component(ll_g) = 'P'
      if (love .gt. 0) component(ll_g) = 'SH'
      
      do k = 1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
      read(15,*)
      read(15,*)(j_con(k), k = 1, jmax)
      read(15,*)(j_wig(k), k = 1, jmax)
      do k = jmin, jmax
         misfit_type(k, ll_g) = j_con(k)
         wavelet_weight(k, ll_g) = real(j_wig(k))
      end do
      do k = jmax+1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
   end do
   close(9)
   close(14)
   close(15)
   
   event_file = 'tele_events.txt'
   inquire(file = event_file, exist = is_file)
   if (is_file) then
      open(12, file=event_file, status='old')
      do ir=1,n_chan
         ll_g = ir+ll_in
         read(12,*)string1, string2, event_sta(ll_g)
      enddo
      close(12)    
   endif
   
   ll_out = ll_in + n_chan
   end subroutine get_body_waves_stations


   subroutine get_surface_waves_stations(ll_in, ll_out)
!
!  Args:
!  ll_in: number of initial channel
!  ll_out: number of final channel
!
   implicit none
   integer ll_in, ll_out, ll_g, io_chan, j_con(11), k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, &
   &  nos(max_stations), io_mod(max_stations), n_chan3, &
   &  io_ns(max_stations), idts, io_ew(max_stations), io_str(max_stations)
   real lat_s(max_stations), lon_s(max_stations), dt, &
   &  ang_ns(max_stations), ang_ew(max_stations), weig(max_stations, 3), j_wig(11), dt_sample, &
   &  dip, rake, theta
   logical, parameter :: cgps=.False., dart=.False.
   character(len=20) filename, string1, string2
   character(len=30) event_file
!   character(len=6) sta_name(max_stations)
   character(len=250) modes
   logical :: is_file

   write(*,*)'Get stations metadata and waveforms for long period surface waves...'
   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
 
   open(9, file='channels_surf.txt', status='old')

   open(15, file='wavelets_surf.txt', status='old')
   read(15,*) jmin, jmax, max_freq
   read(15,'(a)')modes
   read(15,*) n_wave_weight
   idts = 3
 
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
  
   read(9,*) ir_max, n_chan
   read(9,*)
   call error1(ll_in, n_chan)
   call error2(n_wave_weight, n_chan)
!       
!       Here we read the green functions of long period surface waves
!  
   do ir = 1, ir_max    
      ll_g = ir+ll_in
      cgps_channels(ll_g) = .False.
      dart_channels(ll_g) = .False.
      read(9,*) no, sta_name(ll_g), lat_s(ir), lon_s(ir), io_mod(ir), &
      & io_up(ir), io_ns(ir), io_ew(ir), ang_ns(ir), ang_ew(ir), &
      & io_str(ir),(weig(ir, k), k = 1, 3), nos(ir)
      if ((weig(ir, 1) .gt. 0) .and. (io_up(ir) .eq. 1)) n_chan3 = n_chan3 + 1
      if ((weig(ir, 2) .gt. 0) .and. (io_ns(ir) .eq. 1)) n_chan3 = n_chan3 + 1
   end do
   close(9)

   io_chan = 0
   do ir = 1, n_chan
      if (io_up(ir) .eq. 1) then
         io_chan = io_chan + 1
         ll_g = io_chan + ll_in
         weight(ll_g) = weig(ir, 1) / n_chan3
         component(ll_g) = 'R'
      end if
      if (io_ns(ir) .eq. 1) then
         io_chan = io_chan + 1
         ll_g = io_chan + ll_in
         weight(ll_g) = weig(ir, 2) /n_chan3
         component(ll_g) = 'L'
      end if
   end do
   
   filename = 'waveforms_surf.txt'
   filename = trim(filename)
   call get_waveforms(filename, n_chan, ll_in, ll_out) 

   io_chan = 0
   do ir = 1, n_chan
      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      dt_channel(ll_g) = dt_sample
      
      do k = 1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
      read(15,*)
      read(15,*)(j_con(k), k = 1, jmax)
      read(15,*)(j_wig(k), k = 1, jmax)
      do k = jmin, jmax
         misfit_type(k, ll_g) = j_con(k)
         wavelet_weight(k, ll_g) = real(j_wig(k))
      end do
      do k = jmax+1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
   end do
   close(15)   
   
   event_file = 'surf_events.txt'
   inquire(file = event_file, exist = is_file)
   if (is_file) then
      open(12, file=event_file, status='old')
      do ir=1,ir_max
         ll_g = ir+ll_in
         read(12,*)string1, string2, event_sta(ll_g)
      enddo
      close(12)    
   endif
   
   ll_out = ll_in+n_chan
   end subroutine get_surface_waves_stations


   subroutine get_dart_stations(ll_in, ll_out)
!
!  Args:
!  ll_in: number of initial channel
!  ll_out: number of final channel
!
   implicit none
   integer ll_in, ll_out, ll_g, j_con(11), k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, nos, io_mod, n_chan3
   real :: lat_s, lon_s, dt, weig(max_stations), j_wig(11), dt_sample
   logical, parameter :: cgps=.False., dart=.True.
   character(len=20) filename
!   character(len=6) sta_name(max_stations)
!   character(len=3) component(max_stations)

   n_chan3 = 0
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

   read(9,*) ir_max, n_chan
   read(9,*)
   call error1(ll_in, n_chan)
   call error2(n_wave_weight, n_chan)
 
   do ir = 1, ir_max
      ll_g = ir+ll_in
      cgps_channels(ll_g) = .False.
      dart_channels(ll_g) = .True.
      read(9,*) no, sta_name(ll_g), lat_s, lon_s, io_mod, component(ll_g), weig(ir), nos
      if (weig(ir) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   close(9)
   
   filename = 'waveforms_dart.txt'
   filename = trim(filename)
   call get_waveforms(filename, ir_max, ll_in, ll_out) 

   do ir = 1, ir_max

      ll_g = ir+ll_in
      weight(ll_g) = weig(ir) / n_chan3
      dt_channel(ll_g) = dt_sample
      
      do k = 1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
      read(15,*)
      read(15,*)(j_con(k), k = 1, jmax)
      read(15,*)(j_wig(k), k = 1, jmax)
      do k = jmin, jmax
         misfit_type(k, ll_g) = j_con(k)
         wavelet_weight(k, ll_g) = real(j_wig(k))
      end do
      do k = jmax+1, 11
         misfit_type(k, ll_g) = 0
         wavelet_weight(k, ll_g) = 0
      end do
   end do
   close(15)
   ll_out = ll_in + n_chan
   end subroutine get_dart_stations


   subroutine get_waveforms(filename, ir_max, ll_in, ll_out)
!
!  Args:
!  filename: name of file with waveforms
!  ir_max: amount of channels to be read
!  ll_in: number of initial channel
!  ll_out: number of final channel
!
   implicit none
   integer ll_in, ll_out, ll_g, i, ir, ir_max
   real :: dto
   character(len=40) string
   character(len=20) filename
  
   filename = trim(filename)
   open(13, file=filename, status='old')
   do ir = 1, ir_max
      ll_g = ir+ll_in
      t_min(ll_g) = 0
      read(13,*)
      read(13,*)
      read(13,*) string, dto
      read(13,*) string, start(ll_g)
      read(13,*) string, t_max(ll_g)
      read(13,*)
      read(13,*)(obse(i, ll_g), i = 1, start(ll_g))
   enddo
   close(13)
   ll_out = ll_in + ir_max
   end subroutine get_waveforms


   subroutine get_wavelet_obs(wave_obs, wmax, t_max_val)
!
!  Args:
!  wave_obs: name of file with waveforms
!  wmax: amount of channels to be read
!  t_max_val: number of initial channel
!
   implicit none
   real, intent(out) :: wave_obs(wave_pts2, max_stations), wmax(max_stations)
   integer, intent(out) :: t_max_val(max_stations)
   integer nm, i, j, ir, start1, length, n_begin, n_delt
   real cr(wave_pts2), cz(wave_pts2), obser(n_data), &
   &  amp_max, mean
   logical :: cgps, dart
!
! TODO : how to implement this in a more elegant way?
!
   do ir = 1, channels
      start1 = start(ir)
      cgps = cgps_channels(ir)
      dart = dart_channels(ir)
      mean = sum(obse(start1 - 20:start1, ir)) / 20.0
      do i = 1, wave_pts2
         cz(i) = 0.0
         cr(i) = 0.0
         if (i .lt. start1) then
            obser(i) = obse(i, ir)
         else
            obser(i) = 0.0
            if (cgps) obser(i) = mean 
         endif
         
      end do
      do i = 1,nlen
         cr(i) = obser(i)
         cz(i) = 0.0
      enddo

      call wavelet_obs(cr, cz, obser)
      amp_max = 0.0
      nm = 1
      do i = 1, nlen
         if (amp_max .lt. abs(obser(i))) then
            amp_max = abs(obser(i))
            nm = i
         end if
      end do
      if (dart) then
         amp_max = 0.0
         do j = jmin, jmax
            n_begin = 2**(j-1)
            n_delt = nlen/n_begin
            length = int(t_max(ir)/n_delt+0.5)-1
            do i =1, length
               if (amp_max .lt. abs(obser(n_begin+i))) then
                  amp_max = abs(obser(n_begin+i))
                  nm = n_begin+i
               endif
            enddo
         enddo
      endif

      wmax(ir) = amp_max
      t_max_val(ir) = nm
      do i = 1, nlen
         wave_obs(i, ir) = obser(i)/(amp_max)
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
   integer :: j, ir, ir_max, n_begin, n_delt, length
   do ir = 1, channels
      if (weight(ir) .gt. (1e-3)) then
         do j = jmin, jmax
            n_begin = 2**(j-1)
            n_delt = nlen/n_begin
            length = int(t_max(ir)/n_delt+0.5)-1
            if (wavelet_weight(j, ir) .gt. (1e-3)) then
               used_data = used_data + length
            endif
         end do
      endif
   end do
   end subroutine count_wavelets


   subroutine error1(ll_in, n_chan)
!
!  Args:
!  ll_in: number of initial channel
!  n_chan: amount of channels used by specific data type
!
   implicit none
   integer :: ll_in, n_chan
   if (n_chan + ll_in .gt. max_stations) then
      write(*,*)'Error: Maximum allowed number of channels: ', max_stations
      write(*,*)'But amount of channels used for modelling is at least', n_chan + ll_in
      write(*,*)'Please check maximum allowed number of channels'
      stop
   end if
   end subroutine error1
   

   subroutine error2(n_wave_weight, n_chan)
!
!  Args:
!  n_wave_weight: amount of channels specified in wavelet file 
!  n_chan: amount of channels used by specific data type
!
   implicit none
   integer :: n_wave_weight, n_chan
   if (n_wave_weight.ne.n_chan) then
      write(0,*)'n_wave_chan n_chan',n_wave_weight,n_chan
      write(0,*)'Mismatch in amount of stations betweeen wavelet file and stations file'
      stop
   end if
   end subroutine error2


end module get_stations_data      
