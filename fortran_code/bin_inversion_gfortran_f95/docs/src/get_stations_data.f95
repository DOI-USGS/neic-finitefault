module get_stations_data


   use constants, only : nnsta, nnxy, max_seg, npuse, inptd, npth, n_data
   use wavelet_param, only : lnpt, nlen, jmin, jmax, jfmax
   use wavelets, only : meyer_yamada, wavelet_obs
   implicit none
   integer, parameter :: nnsta_tele = 80
   real :: weight(nnsta), wave_obs(npth, nnsta), wmax(nnsta), dt_channel(nnsta)
   integer :: misfit_type(12, nnsta), t_max(nnsta), t_max_val(nnsta)
   real :: wavelet_weight(12, nnsta)


contains


   subroutine get_data(strong, cgnss, body, surf, dart)
   implicit none
   !!
   !!  We load into memory, wavelet transform of observed data, and 
   !!  other properties of stations
   !!
   integer :: ll_in, ll_out
   logical :: strong, cgnss, body, surf, dart
   real erm
   complex z0

   z0 = cmplx(0.0, 0.0)
   erm = 0.0
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
!   select case (io_data)
!      case (1)
!         call get_cgnss_stations(ll_in, ll_out)
!      case (2)
!         call get_strong_motion_stations(ll_in, ll_out)
!      case (3)
!         call get_body_waves_stations(ll_in, ll_out)
!      case (4)
!         call get_strong_motion_stations(ll_in, ll_out)
!         lL_in = ll_out
!         call get_cgnss_stations(ll_in, ll_out)
!      case (5)
!         call get_body_waves_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_cgnss_stations(ll_in, ll_out)
!      case (6)
!         call get_body_waves_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_strong_motion_stations(ll_in, ll_out)
!      case (7)
!         call get_strong_motion_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_surface_waves_stations(ll_in, ll_out)
!      case (8)
!         call get_cgnss_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_surface_waves_stations(ll_in, ll_out)
!      case (9)
!         call get_body_waves_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_surface_waves_stations(ll_in, ll_out)
!      case (10)
!         call get_body_waves_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_strong_motion_stations(ll_in, ll_out)
!         ll_in = lL_out
!         call get_cgnss_stations(ll_in, ll_out)
!      case (11)
!         call get_body_waves_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_surface_waves_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_cgnss_stations(ll_in, ll_out)
!      case (12)
!         call get_body_waves_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_surface_waves_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_strong_motion_stations(ll_in, ll_out)
!      case (13)
!         call get_strong_motion_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_surface_waves_stations(ll_in, ll_out)
!         lL_in = ll_out
!         call get_cgnss_stations(ll_in, ll_out)
!      case (14)
!         call get_body_waves_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_strong_motion_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_cgnss_stations(ll_in, ll_out)
!      case (15)
!         call get_body_waves_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_strong_motion_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_surface_waves_stations(ll_in, ll_out)
!         ll_in = ll_out
!         call get_cgnss_stations(ll_in, ll_out)
!   endselect   
   end subroutine get_data

   
   subroutine get_strong_motion_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, &
   &  nstb, io_chan, nm, j_con(11), i, k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, &
   &  nos(nnsta), io_mod(nnsta), n_chan3
   real cr(inptd), cz(inptd), obser(n_data), &
   &  lat_e, lon_e, lat_s(nnsta), lon_s(nnsta), dt, depth, &
   &  dto, amp_max, weig(nnsta), j_wig(11), df, dt_sample, &
   &  obse(n_data, nnsta), tlen
   character(len=40) string
   character(len=6) sta_name(nnsta)
   character(len=3) component(nnsta)

   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   open(9,file='Readlp.inf',status='old')
   open(13,file='Obser.str',status='unknown')
   open(15,file='Wave.str')
   read(15,*) jmin, jmax, jfmax
   read(15,*)
   read(15,*) n_wave_weight

   read(9,*)
   read(9,*) lat_e, lon_e, depth
   read(9,*) lnpt, dt_sample
   read(9,*)
   nstb = 0
   nlen = 2**lnpt
   dt = dt_sample
   tlen = dt*nlen
   df = 1.0/tlen

   read(9,*) ir_max, n_chan
   read(9,*)
   if (n_wave_weight.ne.n_chan) then
      write(*,*)'n_wave_chan n_chan',n_wave_weight,n_chan
      write(*,*)'Mismatch in amount of stations betweeen wavelet file and stations file'
      stop
   end if
   ir = 0
   io_chan = 0
   df = 1. / ((2 ** lnpt) * dt)
 
   do ir = 1, ir_max
      read(9,*) no, sta_name(ir), lat_s(ir), lon_s(ir), io_mod(ir), component(ir), weig(ir), nos(ir)
      if (weig(ir) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   close(9)
   write(*,*)'n_chan: ', n_chan3

   do ir = 1, ir_max

      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      weight(ll_g) = weig(ir) / n_chan3
      dt_channel(ll_g) = dt_sample
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
   close(13)
   close(15)  
   ll_out = ll_in + n_chan
   end subroutine get_strong_motion_stations

   
   subroutine get_cgnss_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, &
   &  nstb, io_chan, nm, j_con(11), i, k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, &
   &  nos(nnsta), io_mod(nnsta), n_chan3
   real cr(inptd), cz(inptd), obser(n_data), &
   &  lat_e, lon_e, lat_s(nnsta), lon_s(nnsta), dt, depth, &
   &  dto, amp_max, weig(nnsta), j_wig(11), df, dt_sample, &
   &  obse(n_data, nnsta), tlen
   character(len=40) string
   character(len=6) sta_name(nnsta)
   character(len=3) component(nnsta)

   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   open(9,file='Readlp.cgnss',status='old')
   open(13,file='Obser.cgnss',status='unknown')
   open(15,file='Wave.cgnss')
   read(15,*) jmin, jmax, jfmax
   read(15,*)
   read(15,*) n_wave_weight

   read(9,*)
   read(9,*) lat_e, lon_e, depth
   read(9,*) lnpt, dt_sample
   read(9,*)
   nstb = 0
   nlen = 2**lnpt
   dt = dt_sample
   tlen = dt*nlen
   df = 1.0/tlen

   read(9,*) ir_max, n_chan
   read(9,*)
   if (n_wave_weight.ne.n_chan) then
      write(*,*)'n_wave_chan n_chan',n_wave_weight,n_chan
      write(*,*)'Mismatch in amount of stations betweeen wavelet file and stations file'
      stop
   end if
   ir = 0
   io_chan = 0
   df = 1. / ((2 ** lnpt) * dt)
 
   do ir = 1, ir_max
      read(9,*) no, sta_name(ir), lat_s(ir), lon_s(ir), io_mod(ir), component(ir), weig(ir), nos(ir)
      if (weig(ir) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   close(9)

   do ir = 1, ir_max

      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      weight(ll_g) = weig(ir) / n_chan3
      dt_channel(ll_g) = dt_sample
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
   close(13)
   close(15)
   ll_out = ll_in + n_chan
   end subroutine get_cgnss_stations

   
   subroutine get_body_waves_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, &
   &  nstb, nm, j_con(11), i, k, ir, no, n_wave_weight, n_chan, &
   &  nos(nnsta_tele), n_chan3, idts, &
   &  nstaon, mmm(nnsta_tele), love, llove(nnsta_tele), &
   &  idata(nnsta_tele), int1
   real cr(inptd), cz(inptd), obser(n_data), &
   &  lat_sta(nnsta_tele), lon_sta(nnsta_tele), &
   &  dto, amp_max, j_wig(11), dt, &
   &  obse(n_data, nnsta), & 
   &  rang(nnsta_tele), az(nnsta_tele), earth_angle(nnsta_tele), &
   &  disp_or_vel(nnsta_tele), float1, float2 
   character(len=40) string
   character(len=6)  earth(nnsta_tele), sttyp(nnsta_tele), stname(nnsta_tele) 
   character(len=14) fname(nnsta_tele)

   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   open(9,file='Readlp.das',status='old')
   open(13,file='Obser.tele',status='unknown')
   open(14,file='Weight')
   if (jfmax .gt. npuse) then
      write(*,*)'You should stop and check dimension sww,sws'
      stop
   end if
   read(9,*)
   read(9,*) idts, lnpt, dt
   read(9,*)
   read(9,*) nstaon
   nstb = 0
   nlen = 2 ** lnpt

   open(15,file='Wave.tele')
   read(15,*) jmin, jmax, jfmax
   read(15,*) n_wave_weight
   do ir = 1, nstaon
      ll_g = ll_in+ir
      read(14,*) weight(ll_g)
      if (weight(ll_g) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   rewind 14
   write(*,*)'n_chan: ', n_chan3
   n_chan = nstaon

   do ir = 1, nstaon
      read(9,*) nos(ir), earth(ir), sttyp(ir), stname(ir), fname(ir), &
      & rang(ir), az(ir), lat_sta(ir), lon_sta(ir), earth_angle(ir), float1, &
      & mmm(ir), float2, disp_or_vel(ir), llove(ir), int1, idata(ir)
      ll_g = ir+ll_in
      if (idata(ir) .gt. 0 .or. mmm(ir) .eq. 3) cycle
      nstb = nstb+1
      dt_channel(ll_g) = dt
      read(14,*) weight(ll_g)
      weight(ll_g) = weight(ll_g) / n_chan3
      love = llove(ir)
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
   close(13)
   close(14)
   close(15)
   ll_out = ll_in + n_chan
   end subroutine get_body_waves_stations


   subroutine get_surface_waves_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, &
   &  nstb, io_chan, nm, j_con(11), i, k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, &
   &  nos(nnsta), io_mod(nnsta), n_chan3, io_up(nnsta), &
   &  io_ns(nnsta), idts, io_ew(nnsta), io_str(nnsta)
   real cr(inptd), cz(inptd), obser(n_data), &
   &  lat_s(nnsta), lon_s(nnsta), dt, &
   &  ang_ns(nnsta), ang_ew(nnsta), dto, amp_max, &
   &  weig(nnsta, 3), j_wig(11), df, dt_sample, &
   &  obse(n_data, nnsta), dip, rake, theta, tlen
   character(len=40) string
   character(len=6) sta_name(nnsta)
   character(len=250) modes

   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
 
   open(9,file='Readlp.inf_low',status='old')
   open(13,file='Obser.str_low',status='unknown')

   open(15,file='Wave.str_low')
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
   nstb = 0
   nlen = 2**lnpt
   dt = dt_sample
   tlen = dt*nlen
   df = 1.0/tlen
  
   read(9,*) ir_max, n_chan
   read(9,*)
   if (n_wave_weight.ne.n_chan) then
      write(*,*)'n_wave_chan n_chan',n_wave_weight,n_chan
      write(*,*)'Mismatch in amount of stations between wavelet file and stations file'
      stop
   end if
!       
!       Here we read the green functions of long period surface waves
!  
   do ir = 1, ir_max    
      read(9,*) no, sta_name(ir), lat_s(ir), lon_s(ir), io_mod(ir), &
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

   io_chan = 0
   do ir = 1, n_chan
      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      dt_channel(ll_g) = dt_sample
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
   close(13)
   close(15)   
   ll_out = ll_out+n_chan
   end subroutine get_surface_waves_stations


   subroutine get_dart_stations(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, &
   &  nstb, io_chan, nm, j_con(11), i, k, ir, no, &
   &  n_wave_weight, ir_max, n_chan, &
   &  nos(nnsta), io_mod(nnsta), n_chan3
   real cr(inptd), cz(inptd), obser(n_data), &
   &  lat_e, lon_e, lat_s(nnsta), lon_s(nnsta), dt, depth, &
   &  dto, amp_max, weig(nnsta), j_wig(11), df, dt_sample, &
   &  obse(n_data, nnsta), tlen
   character(len=40) string
   character(len=6) sta_name(nnsta)

   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   open(9,file='Readlp.dart',status='old')
   open(13,file='Obser.dart',status='unknown')
   open(15,file='Wave.dart')
   read(15,*) jmin, jmax, jfmax
   read(15,*)
   read(15,*) n_wave_weight

   read(9,*)
   read(9,*) lat_e, lon_e, depth
   read(9,*) lnpt, dt_sample
   read(9,*)
   nstb = 0
   nlen = 2**lnpt
   dt = dt_sample
   tlen = dt*nlen
   df = 1.0/tlen

   read(9,*) ir_max, n_chan
   read(9,*)
   if (n_wave_weight.ne.n_chan) then
      write(*,*)'n_wave_chan n_chan',n_wave_weight,n_chan
      write(*,*)'Mismatch in amount of stations betweeen wavelet file and stations file'
      stop
   end if
   ir = 0
   io_chan = 0
   df = 1. / ((2 ** lnpt) * dt)
 
   do ir = 1, ir_max
      read(9,*) no, sta_name(ir), lat_s(ir), lon_s(ir), io_mod(ir), weig(ir), nos(ir)
      if (weig(ir) .gt. 0) n_chan3 = n_chan3 + 1
   end do
   close(9)

   do ir = 1, ir_max

      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      weight(ll_g) = weig(ir) / n_chan3
      dt_channel(ll_g) = dt_sample
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
   close(13)
   close(15)
   ll_out = ll_in + n_chan
   end subroutine get_dart_stations


end module get_stations_data      
