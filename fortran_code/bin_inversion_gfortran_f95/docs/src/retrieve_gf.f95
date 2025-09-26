!
! TODO: change to velocity if io_v_d = 1
!
module retrieve_gf


   use constants, only : npuse, nnsta, nnxy_m, nnxy, nt, nny, ndis, inptd, max_seg, &
           &  block_far, ltde, n_data, nnpxy, nnpy, nnys, block_stg, pi, twopi, npth
   use model_parameters, only : point_sources, n_seg, v_ref, nxs_sub, nys_sub, delay_seg, &
           &  dxs, dys, nx_p, ny_p, dip_seg, stk_seg, cniu
   use wavelets, only : fourier_coefs, meyer_yamada
   use retrieve_surf_gf, only : get_surf_gf_data, interp_gf, get_surf_gf, npt_bank, dt_bank
   use rad_pattern, only : rad_coef
   use geodesics, only : distaz
   use rise_time, only : fft
   implicit none
   integer, parameter :: nnsta_tele = 80
   complex, allocatable :: green_dip(:, :, :), green_stk(:, :, :)


contains


   subroutine get_gf(strong, cgnss, body, surf, dart)
   !!
   !!  Here, we load into memory the green functions for each subfault, for every used station
   !!
   implicit none
!
   integer ll_in, ll_out
   logical :: strong, cgnss, body, surf, dart
   allocate(green_dip(npuse, nnsta, nnxy_m))
   allocate(green_stk(npuse, nnsta, nnxy_m))

   call fourier_coefs()
   call meyer_yamada()
   ll_in = 0
   ll_out = 0
   if (strong) then
      call get_strong_motion_gf(ll_in, ll_out)
      ll_in = ll_out
   end if
   if (cgnss) then
      call get_cgnss_gf(ll_in, ll_out)
      ll_in = ll_out
   end if
   if (body) then
      call get_body_waves_gf(ll_in, ll_out)
      ll_in = ll_out
   end if
   if (surf) then
      call get_surface_waves_gf(ll_in, ll_out)
      ll_in = ll_out
   end if
   if (dart) then
      call get_dart_gf(ll_in, ll_out)
      ll_in = ll_out
   end if
!   select case (io_data)
!      case (1)
!         call get_cgnss_gf(ll_in, ll_out)
!      case (2)
!         call get_strong_motion_gf(ll_in, ll_out)
!      case (3)
!         call get_body_waves_gf(ll_in, ll_out)
!      case (4)
!         call get_strong_motion_gf(ll_in, ll_out)
!         lL_in = ll_out
!         call get_cgnss_gf(ll_in, ll_out)
!      case (5)
!         call get_body_waves_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_cgnss_gf(ll_in, ll_out)
!      case (6)
!         call get_body_waves_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_strong_motion_gf(ll_in, ll_out)
!      case (7)
!         call get_strong_motion_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_surface_waves_gf(ll_in, ll_out)
!      case (8)
!         call get_cgnss_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_surface_waves_gf(ll_in, ll_out)
!      case (9)
!         call get_body_waves_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_surface_waves_gf(ll_in, ll_out)
!      case (10)
!         call get_body_waves_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_strong_motion_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_cgnss_gf(ll_in, ll_out)
!      case (11)
!         call get_body_waves_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_surface_waves_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_cgnss_gf(ll_in, ll_out)
!      case (12)
!         call get_body_waves_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_surface_waves_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_strong_motion_gf(ll_in, ll_out)
!      case (13)
!         call get_strong_motion_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_surface_waves_gf(ll_in, ll_out)
!         lL_in = ll_out
!         call get_cgnss_gf(ll_in, ll_out)
!      case (14)
!         call get_body_waves_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_strong_motion_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_cgnss_gf(ll_in, ll_out)
!      case (15)
!         call get_body_waves_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_strong_motion_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_surface_waves_gf(ll_in, ll_out)
!         ll_in = ll_out
!         call get_cgnss_gf(ll_in, ll_out)
!   endselect   
   end subroutine get_gf

   
   subroutine get_strong_motion_gf(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, io_v_d, ll_g, ll, &
   &  nstb, io_chan, i, i_s, ir, no, &
   &  n_wave_weight, ir_max, n_chan, &
   &  nos(nnsta), io_mod(nnsta), ixs, iys, n_chan3, lnpt, nlen, &
   &  jmin, jmax, jfmax
   real omega, block, lat_e, lon_e, lat_s(nnsta), lon_s(nnsta), dt, depth, &
   &  weig, df, dt_sample, w, tlen
   complex z0, z
   character(len=80) path
   character(len=6) sta_name(nnsta)
   character(len=80) filename
   character(len=3) component(nnsta), comp
   
   z0 = cmplx(0.0, 0.0)
   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   block = dxs * dys * (1.e-10) 
   write(*,*)'block: ', block 
   write(*,*)'area=',dxs * dys

   open(9,file='Readlp.inf',status='old')
   open(13,file='Obser.str',status='unknown')
   open(15,file='Wave.str')
   read(15,*) jmin, jmax, jfmax
   read(15,'(a)')path
   read(15,*) n_wave_weight

   read(9,*)
   read(9,*) lat_e, lon_e, depth
   read(9,*) lnpt, dt_sample
   read(9,*) io_v_d
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
!       
!       Here we read the green functions of strong motion waves
!
   df = 1. / ((2.0 ** lnpt) * dt)
   nlen = 2 ** lnpt
   tlen = nlen * dt

   do ir = 1, ir_max

      read(9,*) no, sta_name(ir), lat_s(ir), lon_s(ir), io_mod(ir), component(ir), weig, nos(ir)
      comp = component(ir)

      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      if (comp .eq.'HNZ') filename = trim(sta_name(ir))//'1'
      if (comp .eq.'HNN') filename = trim(sta_name(ir))//'2'
      if (comp .eq.'HNE') filename = trim(sta_name(ir))//'3'
      filename = trim(filename)

      open(12,file=filename,status='unknown',access='direct',recl=block_stg)
      ll = 0
!
!       Here, we read the green functions and derivate them
!
      ll_g = ll_in+io_chan
      do i_s = 1, n_seg
         do iys = 1, nys_sub(i_s)
            do ixs = 1, nxs_sub(i_s)
               ll = ll+1
               read(12, rec = ll)(green_dip(i, ll_g, ll), i = 1, jfmax),(green_stk(i, ll_g, ll), i = 1, jfmax)
               do i = 1, jfmax
!
! we eventually shift synthetics in time, in case the fault plane used has a delay
!
                  omega = twopi*(i-1)*df
                  w = -omega*delay_seg(i_s)
                  z = cmplx(0.0, w)
                  z = cexp(z)
                  green_dip(i, ll_g, ll) = green_dip(i, ll_g, ll)*z*block
                  green_stk(i, ll_g, ll) = green_stk(i, ll_g, ll)*z*block
               end do
            end do
         end do
      end do
      close(12)
   end do      
   close(15)
   close(9)
   close(7)
   close(12)
   close(13)
   ll_out = ll_out+n_chan
   end subroutine get_strong_motion_gf

   
   subroutine get_cgnss_gf(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, ll, n_chan3, &
   &  nstb, io_chan, i, i_s, ir, no, lnpt, nlen, &
   &  n_wave_weight, ir_max, n_chan, nos(nnsta), io_mod(nnsta), ixs, iys, &
   &  jmin, jmax, jfmax
   real omega, block, lat_e, lon_e, lat_s(nnsta), lon_s(nnsta), dt, depth, &
   &  weig, df, dt_sample, w, low_freq, tlen
   complex z0, z
   character(len=80) path
   character(len=6) sta_name(nnsta)
   character(len=80) filename
   character(len=3) component(nnsta), comp

   z0 = cmplx(0.0, 0.0)
   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gfs is for Mo = 1e+20
!
   block = dxs * dys * (1.e-10) 
   write(*,*)'block: ', block 
   write(*,*)'area=',dxs * dys

   open(9,file='Readlp.cgnss',status='old')
   open(13,file='Obser.cgnss',status='unknown')
   open(15,file='Wave.cgnss')
   read(15,*) jmin, jmax, jfmax
   read(15,'(a)')path
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
      write(*,*)'Mismatch in number of stations between wavelet file and stations file'
      stop
   end if
   ir = 0
   io_chan = 0
!       
!       Here we read the green functions of strong motion waves
!
   df = 1. / ((2.0 ** lnpt) * dt)
   nlen = 2 ** lnpt
   tlen = nlen * dt
   low_freq = 1.0 / 200.0

   do ir = 1, ir_max

      read(9,*) no, sta_name(ir), lat_s(ir), lon_s(ir), io_mod(ir), component(ir), weig, nos(ir)

      comp = component(ir)

      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      if (comp .eq.'LXZ') filename = trim(sta_name(ir))//'.cgnss.1'
      if (comp .eq.'LXN') filename = trim(sta_name(ir))//'.cgnss.2'
      if (comp .eq.'LXE') filename = trim(sta_name(ir))//'.cgnss.3'
      filename = trim(filename)

      open(12,file=filename,status='unknown',access='direct',recl=block_stg)
      ll = 0
!
!       Here, we read the green functions and derivate them
!
      ll_g = ll_in+io_chan
      do i_s = 1, n_seg
         do iys = 1, nys_sub(i_s)
            do ixs = 1, nxs_sub(i_s)
               ll = ll+1
               read(12, rec = ll)(green_dip(i, ll_g, ll), i = 1, jfmax),(green_stk(i, ll_g, ll), i = 1, jfmax)
               do i = 1, jfmax
                  omega = twopi*(i-1)*df
                  w = -omega*delay_seg(i_s) ! we shift synthetics in time by delay_seg(i_s) seconds
                  z = cmplx(0.0, w)
                  z = cexp(z)
                  green_dip(i, ll_g, ll) = green_dip(i, ll_g, ll)*z*block
                  green_stk(i, ll_g, ll) = green_stk(i, ll_g, ll)*z*block
               end do
            end do
         end do
      end do
      close(12)
   end do
   close(15)
   close(9)
   close(7)
   close(12)
   close(13)
   ll_out = ll_out+n_chan
   end subroutine get_cgnss_gf


   subroutine get_body_waves_gf(ll_in, ll_out)
   implicit none
   integer nstaon, ir, ll_g, k, nsta, n_chan, &
   &  love, ll, i_s, iys, iyp, io_seg, iys_c, iy_c, jf, i, nstb, ipy, npxy, &
   &  nkxy, nxs_c, nys_c, nxp_c, nyp_c, ll_s, kxy, ixs, kpxy, ixp, & 
   &  n_wave_weight, iud(nnsta_tele), &
   &  idata(nnsta_tele), nos(nnsta_tele), llove(nnsta_tele), &
   &  mmm(nnsta_tele), ll_in, ll_out, n_chan3, lnpt, &
   &  nlen, jmin, jmax, jfmax
   real dt, df, ddelt, time, block, w, &
   &  earth_angle(nnsta_tele), hcru(nnsta_tele), &
   &  disp_or_vel(nnsta_tele), rang(nnsta_tele), az(nnsta_tele), &
   &  ttvl(nnsta_tele), &
   &  lat_sta(nnsta_tele), lon_sta(nnsta_tele), tlen
   real, allocatable :: tdel(:,:,:,:)
   complex :: z, z0, wsyn(npth)
   complex :: kahan_y1(inptd), kahan_t1(inptd), kahan_c1(inptd)
   complex :: kahan_y2(inptd), kahan_t2(inptd), kahan_c2(inptd)
   complex, allocatable :: green_dip0(:,:,:,:,:)
   complex, allocatable :: green_stk0(:,:,:,:,:)
!
   character(len=5) earth(nnsta_tele)
   character(len=14) fname(nnsta_tele)
   character(len=6) stname(nnsta_tele), sttyp(nnsta_tele)
   character(len=14) fname4, fname6
   allocate(green_dip0(npuse, nnpy, nnys, max_seg, nnsta_tele))
   allocate(green_stk0(npuse, nnpy, nnys, max_seg, nnsta_tele))
   allocate(tdel(nnpxy, nnxy, max_seg, nnsta_tele))
   
   z0 = cmplx(0.0, 0.0)
   n_chan3 = 0

   open(9,file='Readlp.das',status='old')
   open(13,file='Obser.tele',status='unknown')
   if (jfmax .gt. npuse) then
      write(*,*)'You should stop and check dimension sww,sws'
      stop
   end if
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*) nstaON
   nstb = 0
   block = dxs*dys*(1.0e-10)
   block = block*1.e+4
!  
!  Because the data is micrometer, so The amplitude of block should
!  time e+4 Only in this test
!
   open(15,file='Wave.tele')
   read(15,*) jmin, jmax, jfmax
   read(15,*) n_wave_weight

   tlen = nlen*dt
   npxy = nx_p*ny_p
   do ir = 1, nstaon
      ll_g = ll_in+ir

      read(9,*) nos(ir), earth(ir), sttyp(ir), stname(ir), fname(ir), &
      &  rang(ir), az(ir), lat_sta(ir), lon_sta(ir), earth_angle(ir), ttvl(ir), &
      &  mmm(ir), hcru(ir), disp_or_vel(ir), llove(ir), iud(ir), idata(ir)

      if (idata(ir) .gt. 0 .or. mmm(ir) .eq. 3) cycle
      nstb = nstb+1
      nsta = nstaON
      n_chan = nstaon
      love = llove(ir)
!
!  MM = 0 FOR FAR; MM = 1 FOR UPPER; MM = 3 FOR PNL
!  Here we read the green functions of the teleseismic body waves
!  
      if (love .eq. 0) then
         fname4 = trim(stname(ir))//'.GRE'
         fname6 = trim(stname(ir))//'.TDE'
      else
         fname4 = trim(stname(ir))//'SH.GRE'
         fname6 = trim(stname(ir))//'SH.TDE'
      end if
      open(12,file=fname4,status='unknown',access='direct',recl=block_far)
      open(32,file=fname6,status='unknown',access='direct',recl=ltde)
      LL = 0
      do i_s = 1, n_seg
         if (i_s .le. 10) then
            do iys = 1, nys_sub(i_s)
               do IPY = 1, ny_p
                  LL = LL+1
                  read(12, rec = LL) io_seg, iys_c, iy_c, JF, DT, DF, &
     &  (green_dip0(i, ipy, iys, i_s, ir), i = 1, 2*jfmax)
                  LL = LL+1
                  read(12, rec = LL) io_seg, iys_c, iy_c, JF, LNPT, NLEN, &
     &  (green_stk0(i, ipy, iys, i_s, ir), I = 1, 2*JFmax)
!
!       Sanity checks
!
                  if ((io_seg.ne.i_s) .or. (iys_c.ne.iys) .or. (iy_c.ne.ipy)) then
                     write(*,*)'io_seg vs i_s: ', io_seg,i_s
                     write(*,*)'iys_c, iys: ', iys_c,iys
                     write(*,*)'iy_c, ipy: ', iy_c,ipy
                     write(*,*)"Green function is not matched with fault model"
                     stop
                  end if
                  do i = 1, 2*jfmax
                     green_dip0(i, ipy, iys, i_s, ir) = block*green_dip0(i, ipy, iys, i_s, ir)
                     green_stk0(i, ipy, iys, i_s, ir) = block*green_stk0(i, ipy, iys, i_s, ir)
                  end do
               end do
            end do
            nkxy = nxs_sub(i_s)*nys_sub(i_s)
            read(32, rec = i_s) nxs_c, nys_c, nxp_c, nyp_c, &
     &  ((tdel(k, ll_s, i_s, ir), k = 1, npxy), ll_s = 1, nkxy)
            if ((nxs_c.ne.nxs_sub(i_s)) .or. (nys_c.ne.nys_sub(i_s)) &
     &      .or.(nxp_c.ne.nx_p) .or. (nyp_c.ne.ny_p)) then
                write(*,*)'nxs',nxs_c,nxs_sub(i_s)
                write(*,*)'nys',nys_c,nys_sub(i_s)
                write(*,*)'nxp',nxp_c,nx_p
                write(*,*)'nyp',nyp_c,ny_p
                write(*,'(a)')'Mismatch in amount of point sources or subfaults &
                &between the specified in Fault.time, and those used in the &
                &green functions.'
                stop
             end if
          else
             npxy = nx_p*ny_p
             nkxy = nxs_sub(i_s)*nys_sub(i_s)
             do iys = 1, nys_sub(i_s)
                do ipy = 1, ny_p
                   do i = 1, 2*jfmax
                      green_dip0(i, ipy, iys, i_s, ir) = z0
                      green_stk0(i, ipy, iys, i_s, ir) = z0
                  end do
               end do
            end do
            do ll_s = 1, nkxy
               do k = 1, npxy
                  tdel(k, ll_s, i_s, ir) = 0.0
               end do
            end do
         end if
      end do
      close(12)
      close(32)
   end do
   close(9)
   close(13)
   close(15)
   
   do ir = 1, nstaon
      ll_g = ll_in+ir

      do I = 1, npth
         wsyn(i) = z0
      end do
      LL = 0
      ddelt = 0.0
      z = z0
      do i_s = 1, n_seg
         kxy = 0
         do iys = 1, nys_sub(i_s)
            do ixs = 1, nxs_sub(i_s)
               kxy = kxy+1
               ll = ll+1
               do i = 1, 2*jfmax
                  green_dip(i, ll_g, ll) = z0
                  green_stk(i, ll_g, ll) = z0
               end do
               kahan_y1(:) = z0
               kahan_t1(:) = z0
               kahan_c1(:) = z0
               kahan_y2(:) = z0
               kahan_t2(:) = z0
               kahan_c2(:) = z0
               kpxy = 0
               do iyp = 1, ny_p
                  do ixp = 1, nx_p
                     kpxy = kpxy+1
                     time = tdel(kpxy, kxy, i_s, ir)+&
                  &  point_sources(4, ixp, iyp, ixs, iys, i_s)/v_ref+delay_seg(i_s)
                     do i = 1, 2*jfmax
                        W = -twopi*(i-1)*df*time
                        z = cmplx(0.0, w)
                        z = cexp(z)
!                        green_dip(i, ll_g, ll) = green_dip(i, ll_g, ll)+green_dip0(i, iyp, iys, i_s, ir)*Z
!                        green_stk(i, ll_g, ll) = green_stk(i, ll_g, ll)+green_stk0(i, iyp, iys, i_s, ir)*z
                        kahan_y1(i) = green_dip0(i, iyp, iys, i_s, ir)*z - kahan_c1(i)
                        kahan_t1(i) = green_dip(i, ll_g, ll) + kahan_y1(i)
                        kahan_c1(i) = (kahan_t1(i) - green_dip(i, ll_g, ll)) - kahan_y1(i)
                        green_dip(i, ll_g, ll) = kahan_t1(i)
                        kahan_y2(i) = green_stk0(i, iyp, iys, i_s, ir)*z - kahan_c2(i)
                        kahan_t2(i) = green_stk(i, ll_g, ll) + kahan_y2(i)
                        kahan_c2(i) = (kahan_t2(i) - green_stk(i, ll_g, ll)) - kahan_y2(i)
                        green_stk(i, ll_g, ll) = kahan_t2(i)
                     end do
                  end do
               end do
               do i = 1, 2*jfmax
                  if (int(disp_or_vel(ir)+0.1) .eq. 0) then
                     z = cmplx(1.0, 0.0)
                  elseif (int(disp_or_vel(ir)+0.1) .eq. 1) then
                     w = twopi*(i-1)*df
                     z = cmplx(0.0, w)
                  end if
                  green_dip(i, ll_g, ll) = z*green_dip(i, ll_g, ll)/npxy
                  green_stk(i, ll_g, ll) = z*green_stk(i, ll_g, ll)/npxy
               end do
            end do
         end do
      end do
   end do
   ll_out = ll_out+nstaon
   deallocate(green_stk0)
   deallocate(green_dip0)
   deallocate(tdel)
   end subroutine get_body_waves_gf


   subroutine get_surface_waves_gf(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, nf1, nf2, nf3, nf4, no, &
   &  npp, ix, iy, nxy, j, ll_g, ll, nstb, io_chan, i, k, io_mod(nnsta), &
   &  i_s, ir, i_ch, n_wave_weight, ir_max, n_chan, &
   &  iys, ixs, io_up(nnsta), io_ew(nnsta), io_ns(nnsta), &
   &  n_chan3, jmin, jmax, jfmax, lnpt, nlen

   real*8 :: dip, theta, dis, az, baz, rad_c, coef_v(2, 3), coef_r(2, 5)
   real fint(npth), f1, f2, f3, f4, const_c, tsub(nnpxy), h, omega, &
   &  dist_min, dist_max, depth_sub, dep_min, dep_max, lat_p, lon_p, niu, &
   &  lat_sta, lon_sta, df_bank, tlen_bank, &
   &  time, a, block, lat_e, lon_e, lat_s(nnsta), lon_s(nnsta), dt, rake, &
   &  depth, ang_ns(nnsta), ang_ew(nnsta), df, area, dt_sample, tlen
   
   complex :: kahan_y, kahan_t, kahan_c
   complex sour_sub(npth), green_s(inptd, 10), www, wss, z0

   character(len=250) modes
   character(len=100) surf_gf_bank
   character(len=6) sta_name(nnsta)

   z0 = cmplx(0.0, 0.0)
   n_chan3 = 0
!
!  suppose the ni = u3e+11, then then moment of 1cm*1km^2 
!       subfault is 3e+21. The gf is for Mo = 1e+20 dyne.cm
!  the unit of surface wave is mm.
!
   area = dxs*dys
   block = 1000.0*area*(1.e-10)
  
   open(9,file='Readlp.inf_low',status='old')
   open(13,file='Obser.str_low',status='unknown')
 
   open(15,file='Wave.str_low')
   read(15,*) jmin, jmax, jfmax
   read(15,'(a)')modes
   read(15,*) n_wave_weight
 
   read(9,*)
   read(9,*) lat_e, lon_e, depth
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
   do ir = 1, ir_max    
      read(9,*) no, sta_name(ir), lat_s(ir), lon_s(ir), io_mod(ir), &
      &  io_up(ir), io_ns(ir), io_ew(ir), ang_ns(ir), ang_ew(ir)
   end do
!  int
!  by default frequency window: 0.003 0.004 0.007 0.008
!
   f1 = 0.003
   f2 = 0.004
   f3 = 0.006
   f4 = 0.007
   nf1 = int(f1/df)
   nf2 = int(f2/df)
   nf3 = int(f3/df)
   nf4 = int(f4/df)+1
   do k = 1, nlen/2
      fint(k) = 1.0
      if (k .ge. nf4 .or. k .lt. nf1) fint(k) = 0.0
      if (k .ge. nf1 .and. k .lt. nf2 .and. ((nf2-nf1) .gt. 0)) then
         fint(k) = 0.5*(1.0-cos(pi*(k-nf1)/(nf2-nf1)))
      end if
      if (k .gt. nf3 .and. k .lt. nf4 .and. ((nf4-nf3) .gt. 0)) then
         fint(k) = 0.5*(1.0+cos(pi*(k-nf3)/(nf4-nf3)))
      end if
      fint(k) = fint(k)*block
   end do

   io_chan = 0
   npp = nx_p*ny_p
!       
!       Here we read the green functions of long period surface waves
!
   call get_surf_gf_data(modes, surf_gf_bank)
   tlen_bank = npt_bank*dt_bank
   df_bank = 1.0/tlen_bank
   const_c = dt_bank/dt
   ll = 0
   dist_min = 20.0
   dist_max = 120.0
   do i_s = 1, n_seg
      dip = dip_seg(i_s)
      theta = stk_seg(i_s)
      do iys = 1, nys_sub(i_s)
         nxy = (iys-1)*nxs_sub(i_s)+1
         depth_sub = point_sources(3, 1, int(ny_p/2) + 1, 1, iys, i_s)
         dep_min = point_sources(3, 1, 1, 1, iys, i_s) - 1.0
         dep_max = point_sources(3, 1, ny_p, 1, iys, i_s) + 1.0
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
         do ixs = 1, nxs_sub(i_s)
            ll = ll+1
            nxy = (iys-1)*nxs_sub(i_s)+ixs
            lat_p = point_sources(1, int(nx_p/2) + 1, int(ny_p/2) + 1, ixs, iys, i_s)
            lon_p = point_sources(2, int(nx_p/2) + 1, int(ny_p/2) + 1, ixs, iys, i_s)
            niu = cniu(nxy, i_s)
            do i = 1, jfmax
               omega = (i-1)*df_bank*twopi
               sour_sub(i) = cmplx(0.0, 0.0)
               kahan_y = z0
               kahan_t = z0
               kahan_c = z0
               do iy = 1, ny_p
                  do ix = 1, nx_p
                     h = point_sources(4, ix, iy, ixs, iys, i_s)
                     time = h/v_ref
                     tsub(1) = time+delay_seg(i_s)
                     a = -omega*tsub(1)
                     kahan_y = cmplx(cos(a), sin(a)) - kahan_c
                     kahan_t = sour_sub(i) + kahan_y
                     kahan_c = (kahan_t - sour_sub(i)) - kahan_y
                     sour_sub(i) = kahan_t
!                     sour_sub(i) = sour_sub(i)+cmplx(cos(a), sin(a))
                  end do
               end do
               sour_sub(i) = sour_sub(i)*const_c/npp
            end do
            io_chan = 0
            do ir = 1, ir_max
               lat_sta = lat_s(ir)
               lon_sta = lon_s(ir)
               call distaz(lat_sta, lon_sta, lat_p, lon_p, dis, az, baz)
               dis = dis/111.32d0
               green_s = interp_gf(dis, depth_sub, dist_min, dist_max, dep_min, dep_max)
               do i_ch = 1, 3
                  if (i_ch .eq. 1.and.io_up(ir).ne.1) then
                     cycle
                  elseif (i_ch .eq. 1 .and. io_up(ir) .eq. 1) then
                     rad_c = 0.0
                  end if
                  if (i_ch .eq. 2.and.io_ns(ir).ne.1) then
                     cycle
                  elseif (i_ch .eq. 2 .and. io_ns(ir) .eq. 1) then 
                     rad_c = ang_ns(ir)
                  end if
                  if (i_ch .eq. 3.and.io_ew(ir).ne.1) then
                     cycle
                  elseif (i_ch .eq. 3 .and. io_ew(ir) .eq. 1) then
                     rad_c = ang_ew(ir)
                  end if

                  io_chan = io_chan+1
                  ll_g = io_chan+ll_in
                  call rad_coef(dip, theta, az, rad_c, coef_v, coef_r)
                  do i = 1, jfmax
                     www = cmplx(0.0, 0.0)
                     wss = cmplx(0.0, 0.0)
                     if (i_ch .eq. 1) then
                        do j = 1, 3
                           www = www+coef_v(1, j)*green_s(i, j+5)
                           wss = wss+coef_v(2, j)*green_s(i, j+5)
                        end do
                     else
                        do j = 1, 5
                           www = www+coef_r(1, j)*green_s(i, j)
                           wss = wss+coef_r(2, j)*green_s(i, j)
                        end do
                     end if
                     if (i .le. jfmax .and. i_s .le. 10) then
                        green_stk(i, ll_g, ll) = wss*fint(i)*sour_sub(i)*niu
                        green_dip(i, ll_g, ll) = www*fint(i)*sour_sub(i)*niu
                     else
                        green_dip(i, ll_g, ll) = cmplx(0.0, 0.0)
                        green_stk(i, ll_g, ll) = cmplx(0.0, 0.0)
                     end if
                  end do
               end do
            end do
         end do
      end do
   end do

   close(15)
   close(9)
   close(13)
   ll_out = ll_out+n_chan

   return
   end subroutine get_surface_waves_gf


   subroutine get_dart_gf(ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, io_v_d, ll_g, ll, &
   &  nstb, io_chan, i, i_s, ir, no, ix, iy, npp, &
   &  n_wave_weight, ir_max, n_chan, &
   &  nos(nnsta), io_mod(nnsta), ixs, iys, n_chan3, lnpt, nlen, &
   &  jmin, jmax, jfmax
   real omega, block, lat_e, lon_e, lat_s(nnsta), lon_s(nnsta), dt, depth, &
   &  weig, df, dt_sample, w, tlen, gf_real(inptd), gf_real2(inptd, nnxy_m), &
   &  gf_imag(inptd), h, time, tsub(nnxy_m), a
   complex z0, z, sour_sub(npth)
   complex :: kahan_y, kahan_t, kahan_c
   character(len=80) path
   character(len=6) sta_name(nnsta)
   character(len=80) filename
   
   z0 = cmplx(0.0, 0.0)
   n_chan3 = 0
   npp = nx_p*ny_p
!
   block = (1.e-2)
   write(*,*)'block: ', block 
   write(*,*)'area=',dxs * dys

   open(9,file='Readlp.dart',status='old')
   open(13,file='Obser.dart',status='unknown')
   open(15,file='Wave.dart')
   read(15,*) jmin, jmax, jfmax
   read(15,'(a)')path
   read(15,*) n_wave_weight

   read(9,*)
   read(9,*) lat_e, lon_e, depth
   read(9,*) lnpt, dt_sample
   read(9,*) io_v_d
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
!       
!       Here we read the green functions of strong motion waves
!
   df = 1. / ((2.0 ** lnpt) * dt)
   nlen = 2 ** lnpt
   tlen = nlen * dt
   gf_real2(:, :) = 0.0

   do ir = 1, ir_max

      read(9,*) no, sta_name(ir), lat_s(ir), lon_s(ir), io_mod(ir), weig, nos(ir)

      io_chan = io_chan+1
      ll_g = io_chan+ll_in
      filename = trim(sta_name(ir))//'_gf.txt'
      filename = trim(filename)

      open(12,file=filename)
      do i=1, 181
         read(12,*)(gf_real2(i+10, ll), ll=1,nxs_sub(1)*nys_sub(1))
      enddo
      close(12)
      
      ll = 0
      ix = int(nx_p / 2.0 + 0.51)
      iy = int(ny_p / 2.0 + 0.51)
!
!       Here, we read the green functions
!
      ll_g = ll_in+io_chan
      do i_s = 1, n_seg
         do iys = 1, nys_sub(i_s)
            do ixs = 1, nxs_sub(i_s)
!
!  effect of point sources of each subfault
!
               do i = 1, jfmax
                  omega = (i-1)*df*twopi
                  time = point_sources(5, ix, iy, ixs, iys, i_s)
                  tsub(1) = time + delay_seg(i_s)
                  a = -omega*tsub(1)
!                  sour_sub(i) = cmplx(0.0, 0.0)
                  sour_sub(i) = cmplx(cos(a), sin(a))
!                  kahan_y = z0
!                  kahan_t = z0
!                  kahan_c = z0
!                  do iy = 1, ny_p
!                     do ix = 1, nx_p
!                        h = point_sources(4, ix, iy, ixs, iys, i_s)
!                        time = h/v_ref
!                        tsub(1) = time+delay_seg(i_s)
!                        a = -omega*tsub(1)
!                        kahan_y = cmplx(cos(a), sin(a)) - kahan_c
!                        kahan_t = sour_sub(i) + kahan_y
!                        kahan_c = (kahan_t - sour_sub(i)) - kahan_y
!                        sour_sub(i) = kahan_t
!!                        sour_sub(i) = sour_sub(i)+cmplx(cos(a), sin(a))
!                     end do
!                  end do
!                  sour_sub(i) = sour_sub(i) / npp
               end do
               gf_real(:) = 0.0
               gf_imag(:) = 0.0
               ll = ll+1
               green_stk(:, ll_g, ll) = 0.0
               gf_real(:) = gf_real2(:, ll)
               call fft(gf_real, gf_imag, lnpt, -1.0)
               do i = 1, jfmax
                  green_dip(i, ll_g, ll) = cmplx(gf_real(i), gf_imag(i))
!
! we eventually shift synthetics in time, in case the fault plane used has a delay
!
                  omega = twopi*(i-1)*df
                  w = -omega*(delay_seg(i_s) - 10.0*dt)!9.0*dt - 13)
                  z = cmplx(cos(w), sin(w))
                  green_dip(i, ll_g, ll) = green_dip(i, ll_g, ll)*z*block*sour_sub(i)
               end do
            end do
         end do
      end do
   end do      
   close(15)
   close(9)
   close(7)
   close(13)
   ll_out = ll_out+n_chan
   end subroutine get_dart_gf


   subroutine deallocate_gf()
   deallocate(green_stk)
   deallocate(green_dip)
   end subroutine deallocate_gf


end module retrieve_gf
