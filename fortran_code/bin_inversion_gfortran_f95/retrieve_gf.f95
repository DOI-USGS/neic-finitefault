!
! TODO: change to velocity if io_v_d = 1
!
module retrieve_gf


   use constants, only : npuse, max_stations, max_subfaults, max_subf, nt, nny, ndis, &
   &  wave_pts2, max_seg, block_far, ltde, n_data, max_psources, max_dip_psources, &
   &  max_dip_subfaults, block_stg, pi, twopi, wave_pts
   use wavelet_param, only : get_data_param
   use model_parameters, only : point_sources
   use retrieve_surf_gf, only : get_surf_gf_data, interp_gf, get_surf_gf, npt_bank, &
   &  dt_bank, check_bounds
   use rad_pattern, only : rad_coef
   use geodesics, only : distaz
   use get_stations_data, only : get_properties, disp_or_vel, &
          & idata, mmm, llove, get_event_sta
   implicit none
   integer, parameter :: nnsta_tele = 150
   complex, allocatable :: green_dip(:, :, :), green_stk(:, :, :)
!   complex, allocatable :: green_dip2(:, :, :), green_stk2(:, :, :)
   integer :: segments, nxs_sub(max_seg), nys_sub(max_seg), nx_p, ny_p
   integer :: cum_subfaults(max_seg)
   real :: dip(max_seg), strike(max_seg), delay_seg(max_seg)
   real :: shear(max_subfaults)
   real :: dxs, dys, v_ref, dt_channel(max_stations)
   integer :: lnpt, nlen, jmin, jmax, max_freq, channels, event_sta(max_stations)
   logical :: segment_in_event(max_seg, 10)
   logical :: subfault_in_event(max_subfaults, 10)
   character(len=15) :: sta_name(max_stations)
   character(len=3) :: component(max_stations)
   

contains


   subroutine retrievegf_set_data_properties()
   implicit none
   call get_properties(sta_name, component, dt_channel, channels)
   call get_data_param(lnpt, jmin, jmax, nlen, max_freq)
   call get_event_sta(event_sta)
   end subroutine retrievegf_set_data_properties


   subroutine retrievegf_set_fault_parameters()
   use model_parameters, only : get_shear, get_segments, &
   &  get_subfaults, get_events_segments
   implicit none
   integer :: subfaults
   real :: v_min, v_max
   call get_shear(shear)
   call get_segments(nxs_sub, nys_sub, dip, strike, delay_seg, segments, subfaults, cum_subfaults)
   call get_subfaults(dxs, dys, nx_p, ny_p, v_min, v_max, v_ref)
   call get_events_segments(segment_in_event, subfault_in_event)
   end subroutine retrievegf_set_fault_parameters


   subroutine get_gf(strong, cgps, body, surf, dart, many_events)
!
!  Args:
!  strong: True if strong motion data are used, False otherwise
!  cgps: True if cGPS data are used, False otherwise
!  body: True if body wave data are used, False otherwise
!  surf: True if surface wave data are used, False otherwise
!  dart: True if DART data are used, False otherwise
!  many_events: True if more than one event to be modelled, False otherwise
!
   implicit none
!
   integer first, last
   logical :: strong, cgps, body, surf, dart, many_events
   allocate(green_dip(npuse, max_stations, max_subfaults))
   allocate(green_stk(npuse, max_stations, max_subfaults))

   write(*,*)'Store GF in memory...'
!
!  Here, we load into memory the green functions for each subfault, for every used station
!  
   first = 0
   last = 0
   if (strong) then
      call get_near_field_gf(first, last, many_events, strong, .False.)
      first = last
   end if
   if (cgps) then
      call get_near_field_gf(first, last, many_events, .False., cgps)
      first = last
   end if
   if (body) then
      call get_body_waves_gf(first, last, many_events)
      first = last
   end if
   if (surf) then
      call get_surface_waves_gf(first, last, many_events)
      first = last
   end if
   if (dart) then
      call get_dart_gf(first, last)
      first = last
   end if
   end subroutine get_gf


   subroutine get_near_field_gf(first, last, many_events, strong, cgps)
!
!  Args:
!  first: number of initial channel
!  last: number of final channel
!  many_events: True if more than one event to be modelled, False otherwise
!  strong: True if strong motion data are used, False otherwise
!  cgps: True if cGPS data are used, False otherwise
!
   implicit none
   integer :: first, last, channel, subfault, &
   &  io_chan, i, j, segment, stations, n_chan, &
   &  l, k, event
   real :: omega, factor, start, dt, df, dt_sample, tlen
   complex :: z0, z
   character(len=80) filename, filename2
   character(len=3) comp
   character(len=1) channel2
   character(len=2) event2
   logical :: many_events, strong, cgps
   
   if (strong) write(*,*)'Store strong motion GF in memory...'
   if (cgps) write(*,*)'Store cGPS GF in memory...'
   z0 = cmplx(0.0, 0.0)
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   factor = dxs * dys * (1.e-10) 

   filename = 'channels_strong.txt'
   if (cgps) filename = 'channels_cgps.txt'
   filename = trim(filename)
   open(9, file=filename, status='old')

   read(9,*)
   read(9,*)
   read(9,*) lnpt, dt_sample
   read(9,*)
   nlen = 2**lnpt
   dt = dt_sample
   tlen = dt*nlen
   df = 1.0/tlen

   read(9,*) stations, n_chan
   close(9)
   channel = 0
   io_chan = 0
!       
!       Here we read the green functions of strong motion waves
!

   do i = 1, stations

      io_chan = io_chan+1
      channel = io_chan+first
      comp = component(channel)
      channel2 = comp(3:3)
      filename2 = trim(sta_name(channel))//'.'//comp
      filename2 = trim(filename2)

      open(12, file=filename2, status='old', access='direct',recl=block_stg)
      subfault = 0
!
!       Here, we read the green functions and derivate them
!
      if (many_events) event = event_sta(channel)
      do segment = 1, segments
         do k = 1, nys_sub(segment)
            do l = 1, nxs_sub(segment)
               subfault = subfault+1
               read(12, rec = subfault) &
               & (green_dip(j, channel, subfault), j = 1, max_freq),&
               & (green_stk(j, channel, subfault), j = 1, max_freq)
               do j = 1, max_freq
!
! we eventually shift synthetics in time, in case the fault plane used has a delay
!
                  omega = twopi*(j-1)*df
                  start = -omega*delay_seg(segment)
                  z = cmplx(0.0, start)
                  z = cexp(z)
                  green_dip(j, channel, subfault) = green_dip(j, channel, subfault)*z*factor
                  green_stk(j, channel, subfault) = green_stk(j, channel, subfault)*z*factor
               end do
               if ((many_events) .and. (segment_in_event(segment, event) .eqv. .False.)) then
                  green_dip(:max_freq, channel, subfault) = 0.0
                  green_stk(:max_freq, channel, subfault) = 0.0
               endif
            end do
         end do
      enddo
      close(12)
   end do      
   last = first+n_chan
   end subroutine get_near_field_gf


   subroutine get_body_waves_gf(first, last, many_events)
!
!  Args:
!  first: number of initial channel
!  last: number of final channel
!  many_events: True if more than one event to be modelled, False otherwise
!
   implicit none
   integer nstaon, channel0, channel, subfault, subfault1, psource, &
   &  love, record, segment, iys, iyp, i, psources, &
   &  subfaults, int1, int2, int3, int4, ixs, ixp, first, last, event
   real dt, df, time, factor, w, tlen!, &
   real, allocatable :: tdel(:,:,:,:)
   complex :: z, z0
   complex :: kahan_y1(wave_pts2), kahan_t1(wave_pts2), kahan_c1(wave_pts2)
   complex :: kahan_y2(wave_pts2), kahan_t2(wave_pts2), kahan_c2(wave_pts2)
   complex, allocatable :: green_dip0(:,:,:,:,:)
   complex, allocatable :: green_stk0(:,:,:,:,:)
!
   character(len=14) fname4, fname6
   logical :: many_events
   allocate(green_dip0(npuse, max_dip_psources, max_dip_subfaults, max_seg, nnsta_tele))
   allocate(green_stk0(npuse, max_dip_psources, max_dip_subfaults, max_seg, nnsta_tele))
   allocate(tdel(max_psources, max_subf, max_seg, nnsta_tele))
   
   write(*,*)'Store body waves GF in memory...'
   z0 = cmplx(0.0, 0.0)

   open(9, file='channels_body.txt', status='old')
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*) nstaON
   factor = dxs*dys*(1.0e-10)
   factor = factor*1.e+4
!  
!  Because the data is micrometer, so The amplitude of block should
!  time e+4 Only in this test
!

   dt = dt_channel(first + 1)
   tlen = nlen*dt
   psources = nx_p*ny_p
   do channel0 = 1, nstaon
      channel = first+channel0

      if (idata(channel0) .gt. 0 .or. mmm(channel0) .eq. 3) cycle
      love = llove(channel0)
!
!  MM = 0 FOR FAR; MM = 1 FOR UPPER; MM = 3 FOR PNL
!  Here we read the green functions of the teleseismic body waves
!  
      if (love .eq. 0) then
         fname4 = trim(sta_name(channel))//'.GRE'
         fname6 = trim(sta_name(channel))//'.TDE'
      else
         fname4 = trim(sta_name(channel))//'SH.GRE'
         fname6 = trim(sta_name(channel))//'SH.TDE'
      end if
      open(12, file=fname4, status='old', access='direct', recl=block_far)
      open(32, file=fname6, status='old', access='direct', recl=ltde)
      record = 0
      if (many_events) event = event_sta(channel)
      do segment = 1, segments
         do iys = 1, nys_sub(segment)
            do Iyp = 1, ny_p
               record = record+1
               read(12, rec = record) int1, int2, int3, int4, DT, DF, &
     &  (green_dip0(i, iyp, iys, segment, channel0), i = 1, 2*max_freq)
               record = record+1
               read(12, rec = record) int1, int2, int3, int4, LNPT, NLEN, &
     &  (green_stk0(i, iyp, iys, segment, channel0), I = 1, 2*max_freq)
!
!       Sanity checks
!
               if ((int1.ne.segment) .or. (int2.ne.iys) .or. (int3.ne.iyp)) then
                  write(*,*)'io_seg vs segment: ', int1,segment
                  write(*,*)'iys_c, iys: ', int2,iys
                  write(*,*)'iy_c, ipy: ', int3,iyp
                  write(*,*)"Green function is not matched with fault model"
                  stop
               end if
               do i = 1, 2*max_freq
                  green_dip0(i, iyp, iys, segment, channel0) = &
                  & factor*green_dip0(i, iyp, iys, segment, channel0)
                  green_stk0(i, iyp, iys, segment, channel0) = &
                  & factor*green_stk0(i, iyp, iys, segment, channel0)
               end do
            end do
         end do
         subfaults = nxs_sub(segment)*nys_sub(segment)
         read(32, rec = segment) int1, int2, int3, int4, &
     &  ((tdel(i, subfault, segment, channel0), i = 1, psources), subfault = 1, subfaults)
         if ((int1.ne.nxs_sub(segment)) .or. (int2.ne.nys_sub(segment)) &
     &   .or.(int3.ne.nx_p) .or. (int4.ne.ny_p)) then
            write(*,*)'Subfaults in strike direction',int1,nxs_sub(segment)
            write(*,*)'Subfaults in dip direction',int2,nys_sub(segment)
            write(*,*)'Point sources in strike direction',int3,nx_p
            write(*,*)'Point sources in dip direction',int4,ny_p
            write(*,'(a)')'Mismatch in amount of point sources or subfaults &
            &between the specified in Fault.time, and those used in the &
            &green functions.'
            stop
         end if
         if (many_events .and. (segment_in_event(segment, event) .eqv. .False.)) then
            psources = nx_p*ny_p
            subfaults = nxs_sub(segment)*nys_sub(segment)
            do iys = 1, nys_sub(segment)
               do iyp = 1, ny_p
                  do i = 1, 2*max_freq
                     green_dip0(i, iyp, iys, segment, channel0) = z0
                     green_stk0(i, iyp, iys, segment, channel0) = z0
                  end do
               end do
            end do
            do subfault = 1, subfaults
               do psource = 1, psources
                  tdel(psource, subfault, segment, channel0) = 0.0
               end do
            end do
         endif
      end do
      close(12)
      close(32)
   end do
   close(9)
   
   do channel0 = 1, nstaon
      channel = first+channel0

      z = z0
      subfault = 0
      do segment = 1, segments
         subfault1 = 0
         do iys = 1, nys_sub(segment)
            do ixs = 1, nxs_sub(segment)
               subfault1 = subfault1 + 1
               subfault = subfault+1
               do i = 1, 2*max_freq
                  green_dip(i, channel, subfault) = z0
                  green_stk(i, channel, subfault) = z0
               end do
               kahan_y1(:) = z0
               kahan_t1(:) = z0
               kahan_c1(:) = z0
               kahan_y2(:) = z0
               kahan_t2(:) = z0
               kahan_c2(:) = z0
               psource = 0
               do iyp = 1, ny_p
                  do ixp = 1, nx_p
                     psource = psource+1
                     time = point_sources(4, psource, subfault)/v_ref+delay_seg(segment)
                     time = time + tdel(psource, subfault1, segment, channel0)
                     do i = 1, 2*max_freq
                        W = -twopi*(i-1)*df*time
                        z = cmplx(0.0, w)
                        z = cexp(z)
                        kahan_y1(i) = green_dip0(i, iyp, iys, segment, channel0)*z - kahan_c1(i)
                        kahan_t1(i) = green_dip(i, channel, subfault) + kahan_y1(i)
                        kahan_c1(i) = (kahan_t1(i) - green_dip(i, channel, subfault)) - kahan_y1(i)
                        green_dip(i, channel, subfault) = kahan_t1(i)
                        kahan_y2(i) = green_stk0(i, iyp, iys, segment, channel0)*z - kahan_c2(i)
                        kahan_t2(i) = green_stk(i, channel, subfault) + kahan_y2(i)
                        kahan_c2(i) = (kahan_t2(i) - green_stk(i, channel, subfault)) - kahan_y2(i)
                        green_stk(i, channel, subfault) = kahan_t2(i)
                     end do
                  end do
               end do
               do i = 1, 2*max_freq
                  if (disp_or_vel(channel0) .eq. 0) then
                     z = cmplx(1.0, 0.0)
                  elseif (disp_or_vel(channel0) .eq. 1) then
                     w = twopi*(i-1)*df
                     z = cmplx(0.0, w)
                  end if
                  green_dip(i, channel, subfault) = z*green_dip(i, channel, subfault)/psources
                  green_stk(i, channel, subfault) = z*green_stk(i, channel, subfault)/psources
!                  green_dip2(i, ll, ll_g) = green_dip(i, ll_g, ll)
!                  green_stk2(i, ll, ll_g) = green_stk(i, ll_g, ll)
               end do
            end do
         end do
      end do
   end do
   last = first+nstaon
   deallocate(green_stk0)
   deallocate(green_dip0)
   deallocate(tdel)
   end subroutine get_body_waves_gf


   subroutine get_surface_waves_gf(first, last, many_events)
!
!  Args:
!  first: number of initial channel
!  last: number of final channel
!  many_events: True if more than one event to be modelled, False otherwise
!
   implicit none
   integer first, last, freq1, freq2, freq3, freq4, no, subfault, psource, &
   &  psources, ix, iy, segment_subfault, j, k, io_chan, i, io_mod, &
   &  segment, channel, station, stations, n_chan, event, &
   &  iys, ixs, io_up(max_stations), io_ew(max_stations), io_ns(max_stations)

   real*8 :: dip_segment, theta, dis, az, baz, rad_c, coef_v(2, 3), coef_r(2, 5)
   real filter(wave_pts), freq01, freq02, freq03, freq04, const_c, dist0, omega, &
   &  dist_min, dist_max, depth_sub, dep_min, dep_max, lat_p, lon_p, shear_subfault, &
   &  lat_sta, lon_sta, df_bank, tlen_bank, start0, &
   &  time, factor, lat_s(max_stations), lon_s(max_stations), dt, rake, &
   &  ang_ns(max_stations), ang_ew(max_stations), df, area, dt_sample, tlen
   
   complex :: kahan_y, kahan_t, kahan_c
   complex start(wave_pts), green_s(wave_pts2, 10), green_dip0, green_stk0, z0

   character(len=250) modes
   character(len=100) surf_gf_bank
   character(len=6) sta_name1
   logical :: many_events

   write(*,*)'Store long period surface waves GF in memory...'
   z0 = cmplx(0.0, 0.0)
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gf is for Mo = 1e+20 dyne.cm
!  the unit of surface wave is mm.
!
   area = dxs*dys
   factor = 1000.0*area*(1.e-10)

   open(9, file='surf_filter.txt', status='old')
   read(9,*)freq01, freq02, freq03, freq04
   close(9)
  
   open(9, file='channels_surf.txt', status='old')
   open(15, file='wavelets_surf.txt', status='old')
   read(15,*)
   read(15,'(a)')modes
 
   read(9,*)
   read(9,*) 
   read(9,*) dip_segment, rake, theta, lnpt, dt_sample
   read(9,*)
   if (lnpt.ne.10) then
      write(*,*)"please check input LNPT"
   end if
   nlen = 2**lnpt
   dt = dt_channel(first + 1)!dt_sample
   tlen = dt*nlen
   df = 1.0/tlen
  
   read(9,*) stations, n_chan
   read(9,*)
   do station = 1, stations    
      read(9,*) no, sta_name1, lat_s(station), lon_s(station), io_mod, &
      &  io_up(station), io_ns(station), io_ew(station), ang_ns(station), ang_ew(station)
   end do
!  int
!  by default frequency window: 0.003 0.004 0.007 0.008
!
!   f1 = 0.003
!   f2 = 0.004
!   f3 = 0.006
!   f4 = 0.007
   freq1 = int(freq01/df)
   freq2 = int(freq02/df)
   freq3 = int(freq03/df)
   freq4 = int(freq04/df)+1
   do i = 1, nlen/2
      filter(i) = 1.0
      if (i .ge. freq4 .or. i .lt. freq1) filter(i) = 0.0
      if (i .ge. freq1 .and. i .lt. freq2 .and. ((freq2-freq1) .gt. 0)) then
         filter(i) = 0.5*(1.0-cos(pi*(i-freq1)/(freq2-freq1)))
      end if
      if (i .gt. freq3 .and. i .lt. freq4 .and. ((freq4-freq3) .gt. 0)) then
         filter(i) = 0.5*(1.0+cos(pi*(i-freq3)/(freq4-freq3)))
      end if
      filter(i) = filter(i)*factor!/nlen !filter(k) = filter(k)*block
   end do

   io_chan = 0
   psources = nx_p*ny_p
!       
!       Here we read the green functions of long period surface waves
!
   call get_surf_gf_data(modes, surf_gf_bank)
   tlen_bank = npt_bank*dt_bank
   df_bank = 1.0/tlen_bank
   const_c = dt_bank/dt
   subfault = 0
   dist_min = 20.0
   dist_max = 120.0
   do segment = 1, segments
      dip_segment = dip(segment)
      theta = strike(segment)
      subfault = cum_subfaults(segment) + 1
      dep_min = point_sources(3, 1, subfault)
      subfault = cum_subfaults(segment)
      subfault = subfault + (nys_sub(segment)-1)*nxs_sub(segment) + 1
      psource = (ny_p-1)*nx_p + 1
      dep_max = point_sources(3, psource, subfault)
      call check_bounds(30.0, 90.0, dep_min, dep_max)
      do iys = 1, nys_sub(segment)
         segment_subfault = (iys-1)*nxs_sub(segment)+1
         psource = int(ny_p/2)*nx_p + 1
         subfault = cum_subfaults(segment)
         subfault = subfault + (iys-1)*nxs_sub(segment) + 1
         depth_sub = point_sources(3, psource, subfault)
         dep_min = point_sources(3, 1, subfault) - 1.0
         psource = (ny_p-1)*nx_p + 1
         dep_max = point_sources(3, psource, subfault) + 1.0
         if (dep_min .lt. 4) dep_min = 4
!  
!  Sanity check:
!
         if (dep_max .lt. dep_min) dep_max = 4
!
         call get_surf_gf(surf_gf_bank, dist_min, dist_max, dep_min, dep_max)
!
!       Note that if the tlen is same in both gf bank and used in here
!       values at each frequency are same. In this code, we will let dt_bank = 2 sec
!       and npt_bank = 2048, in contrast, dt_sample = 4 sec, npt_bank = 1024 
!
         subfault = subfault - 1
         do ixs = 1, nxs_sub(segment)
            subfault = subfault + 1
            psource = int(ny_p/2)*nx_p + int(nx_p/2) + 1
            lat_p = point_sources(1, psource, subfault)
            lon_p = point_sources(2, psource, subfault)
            shear_subfault = shear(subfault)
            do i = 1, max_freq
               omega = (i-1)*df_bank*twopi
               start(i) = cmplx(0.0, 0.0)
               kahan_y = z0
               kahan_t = z0
               kahan_c = z0
               do psource = 1, ny_p*nx_p
                  dist0 = point_sources(4, psource, subfault)
                  time = dist0/v_ref + delay_seg(segment)
                  start0 = -omega*time
                  kahan_y = cmplx(cos(start0), sin(start0)) - kahan_c
                  kahan_t = start(i) + kahan_y
                  kahan_c = (kahan_t - start(i)) - kahan_y
                  start(i) = kahan_t
!                  sour_sub(i) = sour_sub(i)+cmplx(cos(a), sin(a))
               end do
               start(i) = start(i)*const_c/psources
            end do
            io_chan = 0
            do station = 1, stations
               lat_sta = lat_s(station)
               lon_sta = lon_s(station)
               call distaz(lat_sta, lon_sta, lat_p, lon_p, dis, az, baz)
               dis = dis/111.32
               green_s = interp_gf(dis, depth_sub, dist_min, dist_max, dep_min, dep_max)
               do i = 1, 3
                  if (i .eq. 1.and.io_up(station).ne.1) then
                     cycle
                  elseif (i .eq. 1 .and. io_up(station) .eq. 1) then
                     rad_c = 0.0
                  end if
                  if (i .eq. 2.and.io_ns(station).ne.1) then
                     cycle
                  elseif (i .eq. 2 .and. io_ns(station) .eq. 1) then 
                     rad_c = ang_ns(station)
                  end if
                  if (i .eq. 3.and.io_ew(station).ne.1) then
                     cycle
                  elseif (i .eq. 3 .and. io_ew(station) .eq. 1) then
                     rad_c = ang_ew(station)
                  end if

                  io_chan = io_chan+1
                  channel = io_chan+first
                  if (many_events) event = event_sta(channel)
                  call rad_coef(dip_segment, theta, az, rad_c, coef_v, coef_r)
                  do j = 1, max_freq
                     green_dip0 = cmplx(0.0, 0.0)
                     green_stk0 = cmplx(0.0, 0.0)
                     if (i .eq. 1) then
                        do k = 1, 3
                           green_dip0 = green_dip0+coef_v(1, k)*green_s(j, k+5)
                           green_stk0 = green_stk0+coef_v(2, k)*green_s(j, k+5)
                        end do
                     else
                        do k = 1, 5
                           green_dip0 = green_dip0+coef_r(1, k)*green_s(j, k)
                           green_stk0 = green_stk0+coef_r(2, k)*green_s(j, k)
                        end do
                     end if
                     if (many_events) then
                        if (j .le. max_freq .and. segment_in_event(segment, event)) then
                           green_stk(j, channel, subfault) = green_stk0*filter(j)*start(j)*shear_subfault
                           green_dip(j, channel, subfault) = green_dip0*filter(j)*start(j)*shear_subfault
                        else
                           green_dip(j, channel, subfault) = cmplx(0.0, 0.0)
                           green_stk(j, channel, subfault) = cmplx(0.0, 0.0)
                        end if
                     else
                        if (j .le. max_freq .and. segment .le. 10) then
                           green_stk(j, channel, subfault) = green_stk0*filter(j)*start(j)*shear_subfault
                           green_dip(j, channel, subfault) = green_dip0*filter(j)*start(j)*shear_subfault
                        else
                           green_dip(j, channel, subfault) = cmplx(0.0, 0.0)
                           green_stk(j, channel, subfault) = cmplx(0.0, 0.0)
                        end if
                     endif
                  end do
               end do
            end do
         end do
      end do
   end do

   close(9)
   last = first+n_chan

   return
   end subroutine get_surface_waves_gf


   subroutine get_dart_gf(first, last)
!
!  Args:
!  first: number of initial channel
!  last: number of final channel
!
   implicit none
   integer :: first, last, channel, subfault, psource, &
   &  i, j, segment, channel_max, n_chan, &
   &  k, l, etc
   real :: omega, dt, df, dt_sample, start, tlen, real1, imag1, time
   complex :: z0, z
   character(len=80) filename
   
   write(*,*)'Store DART GF in memory...'
   z0 = cmplx(0.0, 0.0)
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!

   open(9, file='channels_dart.txt', status='old')

   read(9,*)
   read(9,*)
   read(9,*) lnpt, dt_sample
   read(9,*)
   nlen = 2**lnpt
   dt = dt_sample
   tlen = dt*nlen
   df = 1.0/tlen

   read(9,*) channel_max, n_chan
   close(9)
!       
!       Here we read the green functions of strong motion waves
!
   df = 1. / ((2.0 ** lnpt) * dt)
   nlen = 2 ** lnpt
   tlen = nlen * dt
   psource = int(ny_p/2 + 1)*nx_p + int(nx_p/2 + 1)

   do i = 1, channel_max

      channel = first+i
      filename = trim(sta_name(channel))//'_gf.txt'
      open(12, file=filename, status='old')
!
!       Here, we read the green functions and derivate them
!
      subfault = 0
      do segment = 1, segments
         do k = 1, nys_sub(segment)
            do l = 1, nxs_sub(segment)
               subfault = subfault+1
               time = point_sources(4, psource, subfault)/v_ref
               green_dip(:,channel,subfault) = z0
               green_stk(:,channel,subfault) = z0
               read(12, *)etc, max_freq
               do j = 1, max_freq
                  read(12, *)real1, imag1
                  green_dip(j, channel, subfault) = cmplx(real1, imag1)
               enddo
               do j = 1, max_freq
!
! we eventually shift synthetics in time, in case the fault plane used has a delay
!
                  omega = twopi*(j-1)*df
                  start = -omega*(time+delay_seg(segment))
                  z = cmplx(0.0, start)
                  z = cexp(z)
                  green_dip(j, channel, subfault) = green_dip(j, channel, subfault)*z
               end do
            end do
         end do
      end do
      close(12)
   end do      
   last = first+n_chan
   end subroutine get_dart_gf

   
   subroutine deallocate_gf()
   deallocate(green_stk)
   deallocate(green_dip)
   end subroutine deallocate_gf


end module retrieve_gf
