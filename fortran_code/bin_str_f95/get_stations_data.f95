module get_stations_data


   use constants, only : nnsta, nnxy, max_seg, npuse, inptd, npth, n_data
   use wavelet_param, only : lnpt, nlen, jmin, jmax, jfmax
   use wavelets, only : fourier_coefs, meyer_yamada, wavelet_obs
   implicit none
   integer, parameter :: nnsta_tele = 80
   real :: weight(nnsta), wave_obs(npth, nnsta), wmax(nnsta), dt_channel(nnsta)
   integer :: misfit_type(12, nnsta), t_max(nnsta), t_max_val(nnsta)
   real :: wavelet_weight(12, nnsta)
   character(len=6) :: sta_name1(nnsta), sta_name2(nnsta), sta_name3(nnsta), &
         &     sta_name4(nnsta), sta_name5(nnsta)
   character(len=3) :: component1(nnsta), component2(nnsta), component3(nnsta), &
         &     component4(nnsta), component5(nnsta)
   integer :: mmm(nnsta), llove(nnsta), io_up(nnsta), idata(nnsta)


contains


   subroutine get_data(strong, cgnss, body, surf, dart)
   implicit none
!
!  Here, we load into memory, wavelet transform of observed data, and 
!       other properties of stations
!
   integer :: ll_in, ll_out
   logical :: strong, cgnss, body, surf, dart
   real erm
   complex z0

   z0 = cmplx(0.0, 0.0)
   erm = 0.0
   call fourier_coefs()
   call meyer_yamada()
   ll_in = 0
   ll_out = 0
   dt_channel(:) = 0.0
   if (strong) then
      call get_strong_motion_stations(ll_in, ll_out)
      ll_in = ll_out
   end if
   if (cgnss) then
      call get_cgnss_stations(ll_in, ll_out)
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
   end subroutine get_data

   
   subroutine get_strong_motion_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, io_chan, j_con(11), k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, &
   &  nos(nnsta), io_mod(nnsta), n_chan3
   real dt, weig(nnsta), j_wig(11), dt_sample, lat_s, lon_s
   character(len=20) filename
!   character(len=6) sta_name(nnsta)
!   character(len=3) component(nnsta)

   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   open(9,file='Readlp.inf',status='old')
   open(15,file='Wave.str', status='old')
   read(15,*) jmin, jmax, jfmax
   read(15,*)
   read(15,*) n_wave_weight

   read(9,*)
   read(9,*) 
   read(9,*) lnpt, dt_sample
   read(9,*)
   nlen = 2**lnpt
   dt = dt_sample

   read(9,*) ir_max, n_chan
   read(9,*)
   call error1(ll_in, n_chan)
   if (n_wave_weight.ne.n_chan) then
      write(*,*)'n_wave_chan n_chan',n_wave_weight,n_chan
      write(*,*)'Mismatch in amount of stations betweeen wavelet file and stations file'
      stop
   end if
   ir = 0
   io_chan = 0
 
   do ir = 1, ir_max
      read(9,*) no, sta_name1(ir), lat_s, lon_s, io_mod(ir), component1(ir), weig(ir), nos(ir)
      if (weig(ir) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   close(9)
   write(*,*)'n_chan: ', n_chan3

   filename = 'Obser.str'
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
   ll_out = ll_in + n_chan
   end subroutine get_strong_motion_stations

   
   subroutine get_cgnss_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, io_chan, j_con(11), k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, &
   &  nos(nnsta), io_mod(nnsta), n_chan3
   real :: lat_s, lon_s, dt, weig(nnsta), j_wig(11), dt_sample
   character(len=20) filename
!   character(len=6) sta_name(nnsta)
!   character(len=3) component(nnsta)

   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   open(9, file='Readlp.cgnss', status='old')
   open(15, file='Wave.cgnss', status='old')
   read(15,*) jmin, jmax, jfmax
   read(15,*)
   read(15,*) n_wave_weight

   read(9,*)
   read(9,*)
   read(9,*) lnpt, dt_sample
   read(9,*)
   nlen = 2**lnpt
   dt = dt_sample

   read(9,*) ir_max, n_chan
   read(9,*)
   call error1(ll_in, n_chan)
   if (n_wave_weight.ne.n_chan) then
      write(*,*)'n_wave_chan n_chan',n_wave_weight,n_chan
      write(*,*)'Mismatch in amount of stations betweeen wavelet file and stations file'
      stop
   end if
   ir = 0
   io_chan = 0
 
   do ir = 1, ir_max
      read(9,*) no, sta_name2(ir), lat_s, lon_s, io_mod(ir), component2(ir), weig(ir), nos(ir)
      if (weig(ir) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   close(9)

   filename = 'Obser.cgnss'
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
   ll_out = ll_in + n_chan
   end subroutine get_cgnss_stations

   
   subroutine get_body_waves_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, j_con(11), k, ir, n_wave_weight, n_chan, &
   &  nos, n_chan3, idts, nstaon, love, int1
   real lat_sta, lon_sta, j_wig(11), dt, & 
   &  rang, az, earth_angle, disp_or_vel(nnsta), float1, float2
   character(len=20) filename
   character(len=6)  earth, sttyp 
   character(len=14) fname

   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   open(9, file='Readlp.das', status='old')
   open(14, file='Weight', status='old')
   if (jfmax .gt. npuse) then
      write(*,*)'You should stop and check dimension sww,sws'
      stop
   end if
   read(9,*)
   read(9,*) idts, lnpt, dt
   read(9,*)
   read(9,*) nstaon
   call error1(ll_in, nstaon)
   nlen = 2 ** lnpt

   open(15, file='Wave.tele', status='old')
   read(15,*) jmin, jmax, jfmax
   read(15,*)
   read(15,*) n_wave_weight
   do ir = 1, nstaon
      ll_g = ll_in+ir
      read(14,*) weight(ll_g)
      if (weight(ll_g) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   rewind 14
   write(*,*)'n_chan: ', n_chan3
   n_chan = nstaon 
   
   filename = 'Obser.tele'
   filename = trim(filename)
   call get_waveforms(filename, nstaon, ll_in, ll_out) 

   do ir = 1, nstaon
      read(9,*) nos, earth, sttyp, sta_name3(ir), fname, &
      & rang, az, lat_sta, lon_sta, earth_angle, float1, &
      & mmm(ir), float2, disp_or_vel(ir), llove(ir), int1, idata(ir)
      ll_g = ir+ll_in
      if (idata(ir) .gt. 0 .or. mmm(ir) .eq. 3) cycle
      dt_channel(ll_g) = dt
      read(14,*) weight(ll_g)
      weight(ll_g) = weight(ll_g) / n_chan3
      love = llove(ir)
      
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
   ll_out = ll_in + n_chan
   end subroutine get_body_waves_stations


   subroutine get_surface_waves_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, io_chan, j_con(11), k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, &
   &  nos(nnsta), io_mod(nnsta), n_chan3, &
   &  io_ns(nnsta), idts, io_ew(nnsta), io_str(nnsta)
   real lat_s(nnsta), lon_s(nnsta), dt, &
   &  ang_ns(nnsta), ang_ew(nnsta), weig(nnsta, 3), j_wig(11), dt_sample, &
   &  dip, rake, theta
   character(len=20) filename
!   character(len=6) sta_name(nnsta)
   character(len=250) modes

   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
 
   open(9, file='Readlp.inf_low', status='old')

   open(15, file='Wave.str_low', status='old')
   read(15,*) jmin, jmax, jfmax
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
  
   read(9,*) ir_max, n_chan
   read(9,*)
   call error1(ll_in, n_chan)
   if (n_wave_weight.ne.n_chan) then
      write(*,*)'n_wave_chan n_chan',n_wave_weight,n_chan
      write(*,*)'Mismatch in amount of stations between wavelet file and stations file'
      stop
   end if
!       
!       Here we read the green functions of long period surface waves
!  
   do ir = 1, ir_max    
      read(9,*) no, sta_name4(ir), lat_s(ir), lon_s(ir), io_mod(ir), &
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
      end if
      if (io_ns(ir) .eq. 1) then
         io_chan = io_chan + 1
         ll_g = io_chan + ll_in
         weight(ll_g) = weig(ir, 2) /n_chan3
      end if
   end do
   
   filename = 'Obser.str_low'
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
   ll_out = ll_out+n_chan
   end subroutine get_surface_waves_stations


   subroutine get_dart_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, j_con(11), k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, nos, io_mod, n_chan3
   real :: lat_s, lon_s, dt, weig(nnsta), j_wig(11), dt_sample
   character(len=20) filename
!   character(len=6) sta_name(nnsta)
!   character(len=3) component(nnsta)

   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   open(9, file='Readlp.dart', status='old')
   open(15, file='Wave.dart', status='old')
   read(15,*) jmin, jmax, jfmax
   read(15,*)
   read(15,*) n_wave_weight

   read(9,*)
   read(9,*)
   read(9,*) lnpt, dt_sample
   read(9,*)
   nlen = 2**lnpt
   dt = dt_sample

   read(9,*) ir_max, n_chan
   read(9,*)
   call error1(ll_in, n_chan)
   if (n_wave_weight.ne.n_chan) then
      write(*,*)'n_wave_chan n_chan',n_wave_weight,n_chan
      write(*,*)'Mismatch in amount of stations betweeen wavelet file and stations file'
      stop
   end if
 
   do ir = 1, ir_max
      read(9,*) no, sta_name5(ir), lat_s, lon_s, io_mod, component5(ir), weig(ir), nos
      if (weig(ir) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   close(9)
   
   filename = 'Obser.dart'
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
   implicit none
   integer ll_in, ll_out, ll_g, nm, i, ir, no, &
   &  ir_max
   real cr(inptd), cz(inptd), obser(n_data), &
   &  dto, amp_max, obse(n_data, nnsta)
   character(len=40) string
   character(len=20) filename
  
   filename = trim(filename)
   open(13, file=filename, status='old')
   do ir = 1, ir_max
      ll_g = ir+ll_in
      read(13,*)
      read(13,*)
      read(13,*) string, dto
      read(13,*) string, no
      read(13,*) string, t_max(ll_g)
      read(13,*)
      read(13,*)(obse(i, ll_g), i = 1, no)
      do i = 1, inptd
         cz(i) = 0.0
         cr(i) = 0.0
         if (i .lt. no) then
            obser(i) = obse(i, ll_g)
         else
            obser(i) = 0.0
         end if
      end do

      call wavelet_obs(cr, cz, obser)
      amp_max = 0.0
      nm = 1
      do i = 1, nlen
         if (amp_max .lt. abs(obser(i))) then
            amp_max = abs(obser(i))
            nm = i
         end if
      end do

      wmax(ll_g) = amp_max
      t_max_val(ll_g) = nm
      do i = 1, nlen
         wave_obs(i, ll_g) = obser(i)/(amp_max)
      end do
   enddo
   close(13)
   ll_out = ll_in + ir_max
   end subroutine get_waveforms


   subroutine error1(ll_in, n_chan)
   integer :: ll_in, n_chan
   if (n_chan + ll_in .gt. nnsta) then
      write(*,*)'Error: Maximum allowed number of channels: ', nnsta
      write(*,*)'But amount of channels used for modelling is at least', n_chan + ll_in
      write(*,*)'Please check maximum allowed number of channels'
      stop
   end if
   end subroutine error1


end module get_stations_data      
