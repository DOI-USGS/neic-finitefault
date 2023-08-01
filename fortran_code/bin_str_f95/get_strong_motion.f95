!c
!c	For strong motion , we need use geometry coordanation directly.
!c
program get_strong_motion


   use constants, only : dpi, twopi, max_dip_psources, max_stk_psources, max_dip_subfaults, max_stk_subfaults, &
             &     max_stk_psources2, max_dip_psources2, wave_pts2, max_seg, block_stg, nt
   use vel_model_data, only : new_vel_s, new_dens, s_layer, read_vel_model, update_model
   use wave_travel, only : trav
   use retrieve_gf, only : dt_gfs, lnpt_gfs, t_cor, get_gf_data, get_gf, interp_gf, deallocate_gf
   use geodesics, only : distaz 
   use fast_fourier_trans, only : fft, fourier_coefs
   use bpfilter, only : bandpassfilter
   use rad_pattern, only : rad_coef
   implicit none
   character(len=100) :: vel_model, gf_file, gf_bank, filter_file, wave_file, stat_file
   character(len=70) :: string1, input
   character(len=200) :: directory,filterfile,risetimefile,pointsourcefile
   character(len=200) :: channelsfile,waveformsfile,responsefile
   character(len=10) :: sta_name(200)
   character(len=12) :: filename
   character(len=3) :: comp
   character(len=1) :: channel2
   real*8 :: dip, theta, az, baz, ang_d, coef_v(2,3), coef_r(2,5)
   real :: dep_p, dep_min, dep_max
   real :: lat_e, lon_e, dep_e
   real :: fau_mod(7, max_stk_psources, max_dip_psources, max_stk_subfaults, max_dip_subfaults, max_seg)
   real :: lat_sta, lon_sta, lat_p, lon_p
   real*8 :: azi(max_stk_psources2, max_dip_psources2), bazi(max_stk_psources2, max_dip_psources2)
   real*8 :: distances(max_stk_psources2, max_dip_psources2), dis
   real :: time(max_stk_psources2,max_dip_psources2), dist(max_stk_psources2), tmin(max_stk_psources2)
   real :: dt, kahan_y1(wave_pts2), kahan_y2(wave_pts2), kahan_t1(wave_pts2), kahan_t2(wave_pts2), &
   &    kahan_c1(wave_pts2), kahan_c2(wave_pts2)
   real*8 :: df, w
   real :: t1, real_dip(wave_pts2), imag_dip(wave_pts2), niu
   real :: real_stk(wave_pts2), imag_stk(wave_pts2)
   real :: niu_fault(max_dip_psources, max_dip_subfaults, max_seg)
   real :: depth_sources(max_dip_psources, max_dip_subfaults, max_seg)
   real :: low_freq, high_freq
   integer :: i, j, nsta
   integer :: ir_max, n_chan, ir, nx, ky, kx, LL, no
   integer :: iys, ixs, iyp, ixp, nxp, nyp
   integer :: npt, lnpt, k, io_mod(200), nn_comp
   integer :: jfmax
   integer :: n_pos, n_delt, nleft
   real :: green_s(nt, 8), v_ref, green
   real :: ww, ws, green_dip(wave_pts2), green_stk(wave_pts2)
   complex :: ssww(wave_pts2), ssws(wave_pts2)
   complex*16 :: z, z1
   integer :: n_seg, i_seg
   integer :: nxs_sub(max_seg), nys_sub(max_seg)
   real :: float1, float2
   integer :: integer1, integer2

   real :: dip_sub(max_seg), stk_sub(max_seg)
   real :: c_depth
   logical :: disp, is_file

   call getarg(1, input)
   call getarg(2, directory)
   disp = (input.eq.'cgps')
   write(*,'(/A/)')'METHOD TO STORE NEAR-FIELD GF'
!
!	Input fault data
!
   call fourier_coefs()
   write(*,*)'Get fault segments data...'
   risetimefile=trim(directory)//'fault&rise_time.txt'
   open(22, file=risetimefile, status='old')
   read(22, *)integer1, integer2, c_depth
   read(22, *)n_seg, float1, float2, nxp, nyp
   read(22, *)float1, float2, integer1, v_ref
   do i_seg = 1, n_seg
      read(22, *)integer1, dip_sub(i_seg),stk_sub(i_seg)
      read(22, *)nxs_sub(i_seg),nys_sub(i_seg)
      do iys = 1, nys_sub(i_seg)
         do ixs = 1, nxs_sub(i_seg)
            read(22,*)
         enddo
      enddo
   enddo
   close(22)
   pointsourcefile=trim(directory)//'point_sources.txt'
   open(22, file=pointsourcefile, status='old')
   do i_seg=1,n_seg
      read(22,*)
      do iys = 1, nys_sub(i_seg)
         do ixs = 1, nxs_sub(i_seg)
            do iyp = 1, nyp
               do ixp = 1, nxp
                  read(22,*)(fau_mod(k, ixp, iyp, ixs, iys, i_seg),k=1,7)
               enddo
            enddo
         enddo
      enddo
   enddo
   close(22)

   filter_file = 'filtro_strong.txt'
   inquire( file = trim(directory)//'filtro_cgps.txt', exist = is_file )
   if (is_file .and. disp) filter_file = 'filtro_cgps.txt'
   filterfile=trim(directory)//filter_file
   open(1, file=filterfile, status='old')
   read(1,*)string1, low_freq, high_freq
   close(1)
   write(*,*)high_freq
   if (disp) low_freq = 0.0
!
!	Load vel. model into memory
! 
   gf_file = trim(directory)//'Green_strong.txt'
   if (disp) gf_file = trim(directory)//'Green_cgps.txt'
   call get_gf_data(gf_file, vel_model, gf_bank)
   write(*,*)'Get velocity model...'
   call read_vel_model(trim(directory)//'vel_model.txt')
   dep_min = 1.e+10
   dep_max = 0.0
   do i_seg=1,n_seg
      do iys = 1, nys_sub(i_seg)
         do iyp = 1, nyp
            dep_p = fau_mod(3, 1, iyp, 1, iys, i_seg)
            if(dep_p.lt.dep_min)dep_min = dep_p
            if(dep_p.gt.dep_max)dep_max = dep_p
            call update_model(dep_p)
            niu = new_vel_s(s_layer)*new_vel_s(s_layer)*new_dens(s_layer) * (1.e+10)
            niu_fault(iyp, iys, i_seg) = niu
            depth_sources(iyp, iys, i_seg) = dep_p
         enddo
      enddo
   enddo
   write(*,*)'Get GF bank...'
   call get_gf(gf_bank, dep_min, dep_max)
!
   wave_file = 'wavelets_strong.txt'
   if (disp) wave_file = 'wavelets_cgps.txt'
   waveformsfile=trim(directory)//wave_file
   open(13, file=waveformsfile, status='old')
   read(13,*)integer1, integer2, jfmax
   read(13,*)
   read(13,*)nsta
   if(nsta.lt.1)then
      write(*,*)"In wavelet file, you have not specified any stations"
      stop
   endif
   close(13)
!
!	Load hypocenter location
! 
   write(*,*)'Get near field stations data...'
   stat_file = 'channels_strong.txt'
   if (disp) stat_file = 'channels_cgps.txt'
   channelsfile=trim(directory)//stat_file
   open(12, file=channelsfile, status='old')
   read(12,*)
   read(12,*)lat_e, lon_e, dep_e
   read(12,*)lnpt,dt
   read(12,*)
   if(abs(c_depth-dep_e).gt.1e-4)then
      write(*,*)"c_depth,dep_e"
      write(*,*)c_depth,dep_e
      write(*,*)'Mismatch between hypocenter depth in Readlp.inf and in &
     & Fault.time'
      stop
   endif
   npt = 2 ** lnpt
   nleft = npt - 80
   if(npt.gt.nt)then
      write(*,*)"The length of green's function must be shorten"
      stop
   endif
   df = 1.d0 / npt / dt
   if(abs(lnpt_gfs-lnpt).gt.1e-6.or.abs(dt-dt_gfs).gt.1.0e-6)then
      write(0,*)"Please keep the consistency between Readlp.inf and Gfs bank"
      write(0,*)"In Readlp.inf dt= ",dt, " Lnpt= ",lnpt
      write(0,*)"In Gfs        dt= ",dt_gfs," Lnpt= ",lnpt_gfs
      stop
   endif

   read(12,*)Ir_max,n_chan
   read(12,*)
!   write(*,*)'Get near field stations data...'
   write(*,*)'Begin to store the response for each channel...'
!
!	Load GF for each station and channel
!
   do ir = 1, ir_max
      read(12,*)no, sta_name(ir), lat_sta, lon_sta, io_mod(ir), comp
      write(*,*)'Store response for station ', sta_name(ir), ', channel ', comp, '...'
      call distaz(lat_sta,lon_sta,lat_e,lon_e,dis,az,baz)
!      write(*,*)'az=',az,dis,baz,lat_sta,lon_sta,lat_e,lon_e
      channel2 = comp(3:3)
!      write(*,*)scan(comp, 'Z')
      if (channel2 .eq. 'Z') then
         nn_comp = 14
      elseif (channel2 .eq. 'N') then
         nn_comp = 15
      elseif (channel2 .eq. 'E') then
         nn_comp = 16
      endif
      filename = trim(sta_name(ir))//'.'//comp
      open(20,file=trim(directory)//filename,status='unknown',access='direct',recl=block_stg)
      ll = 0
      ky = 0
      do i_seg = 1,n_seg
         dip = dip_sub(i_seg)
         theta = stk_sub(i_seg)
         do iys = 1, nys_sub(i_seg)
!
!	 Store azimuth, back-azimuth, distance, and arrivals, of point sources
!
            do iyp = 1, nyp
               ky = (iys - 1) * nyp + iyp
               kx = 0
               do ixs = 1, nxs_sub(i_seg)
                  do ixp = 1, nxp
                     kx = kx + 1
                     lat_p = fau_mod(1, ixp, iyp, ixs, iys, i_seg)
                     lon_p = fau_mod(2, ixp, iyp, ixs, iys, i_seg)
                     call distaz(lat_sta, lon_sta, lat_p, lon_p, dis, az, baz)
                     if (dis.lt.0.25) dis = 0.25d0
                     dist(kx) = dis
                     distances(kx, ky) = dis
                     azi(kx,ky) = az
                     bazi(kx,ky) = baz
                  enddo
               enddo
               nx = nxs_sub(i_seg) * nxp
               dep_p = depth_sources(iyp, iys, i_seg)
               call update_model(dep_p)
               call trav(dist, nx, tmin)
               time(:nx,ky) = tmin(:nx)
            enddo
!
!	Calculate subfault response
!	First save the file: vertical:  STA.1 - STA.2 - STA.3 
!
            do ixs = 1, nxs_sub(i_seg)
               LL = LL + 1
               green_dip(:) = 0.0
               green_stk(:) = 0.0
               kahan_y1(:) = 0.0
               kahan_y2(:) = 0.0
               kahan_t1(:) = 0.0
               kahan_t2(:) = 0.0
               kahan_c1(:) = 0.0
               kahan_c2(:) = 0.0
               do iyp = 1, nyp
                  niu = niu_fault(iyp, iys, i_seg)
                  dep_p = depth_sources(iyp, iys, i_seg)
                  ky = (iys - 1) * nyp + iyp
                  do ixp = 1, nxp
                     kx = (ixs - 1) * nxp + ixp
                     dis = distances(kx, ky)
                     az = azi(kx,ky)
                     baz = bazi(kx,ky)   
                     t1 = fau_mod(4, ixp, iyp, ixs, iys, i_seg) / v_ref
!                     t1 = min(t1, t_latest)
                     t1 = t1 + time(kx, ky)
                     n_delt = int(t1 / dt +0.5)
                     green_s = interp_gf(dis, dep_p, dep_min, dep_max) 
                     ang_d = (baz - 180.0d0)*dpi
                     call rad_coef(dip,theta,az,ang_d,nn_comp,coef_v,coef_r)
!
! n_delt: position of array for which the wave arrives to the station from the given subfault.
! n_pos = i - n_delt, means we initialize the green function when p-wave from the subfault arrives to the station
!
                     if(nn_comp.eq.14)then
                        do i = n_delt + 1, npt
                           ww = 0.0
                           ws = 0.0
                           n_pos = i - n_delt
                           do j = 1, 3
                              green = green_s(n_pos, j + 5) * niu
                              ww = ww + coef_v(1, j) * green
                              ws = ws + coef_v(2, j) * green
                           enddo
                           kahan_y1(i) = ww - kahan_c1(i)
                           kahan_t1(i) = green_dip(i) + kahan_y1(i)
                           kahan_c1(i) = (kahan_t1(i) - green_dip(i)) - kahan_y1(i)
                           green_dip(i) = kahan_t1(i)
                           kahan_y2(i) = ws - kahan_c2(i)
                           kahan_t2(i) = green_stk(i) + kahan_y2(i)
                           kahan_c2(i) = (kahan_t2(i) - green_stk(i)) - kahan_y2(i)
                           green_stk(i) = kahan_t2(i)
!                           green_dip(i) = green_dip(i) + ww
!                           green_stk(i) = green_stk(i) + ws
                        enddo
                     else
                        do i = n_delt + 1, npt
                           ww = 0.0
                           ws = 0.0
                           n_pos = i - n_delt
                           do j = 1,5
                              green = green_s(n_pos, j) * niu
                              ww = ww + coef_r(1, j) * green
                              ws = ws + coef_r(2, j) * green
                           enddo
                           kahan_y1(i) = ww - kahan_c1(i)
                           kahan_t1(i) = green_dip(i) + kahan_y1(i)
                           kahan_c1(i) = (kahan_t1(i) - green_dip(i)) - kahan_y1(i)
                           green_dip(i) = kahan_t1(i)
                           kahan_y2(i) = ws - kahan_c2(i)
                           kahan_t2(i) = green_stk(i) + kahan_y2(i)
                           kahan_c2(i) = (kahan_t2(i) - green_stk(i)) - kahan_y2(i)
                           green_stk(i) = kahan_t2(i)
!                           green_dip(i) = green_dip(i) + ww
!                           green_stk(i) = green_stk(i) + ws
                        enddo
                     endif
                  enddo
               enddo
!
! synthetic in dip and strike component
!
               imag_dip(:) = 0.0
               real_dip(:) = 0.0
               real_dip(:nleft) = green_dip(:nleft) / (nxp * nyp)

               call bandpassfilter(real_dip, wave_pts2, dt, 4, 2, low_freq, high_freq)
               call fft(real_dip, imag_dip, lnpt, -1)
               imag_stk(:) = 0.0
               real_stk(:) = 0.0
               real_stk(:nleft) = green_stk(:nleft) / (nxp * nyp)
               
               call bandpassfilter(real_stk, wave_pts2, dt, 4, 2, low_freq, high_freq)
               call fft(real_stk, imag_stk, lnpt, -1)
               w = t_cor*twopi*df
               z1 = cmplx(cos(w), sin(w), kind(1d0))
               z = cmplx(1.d0, 0.d0, kind(1d0))
               do i = 1, jfmax
                  ssww(i) = cmplx(real_dip(i), imag_dip(i)) * z
                  ssws(i) = cmplx(real_stk(i), imag_stk(i)) * z
                  z = z*z1
               enddo
               write(20, rec=ll)(ssww(i),i=1,jfmax),(ssws(i),i=1,jfmax)
            enddo
         enddo
      enddo
      close(20)
   enddo
   close(12)
   call deallocate_gf()

   write(*,'(/A/)')'END METHOD TO STORE NEAR-FIELD GF'


end program get_strong_motion
