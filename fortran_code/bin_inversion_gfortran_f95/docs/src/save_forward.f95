module save_forward


   use constants, only : max_seg, nnxy, inptd, npth, nnsta, n_data, nnpxy, twopi, &
           &      nnpy, nnys, npuse, pi, mmsou, dpi 
   use model_parameters, only : nxs_sub, nys_sub, ta0, dta, msou, n_seg, dxs, dys, &
           &      nx_p, ny_p
   use wavelet_param, only : lnpt, jfmax, nlen 
   use get_stations_data, only : dt_channel 
   use retrieve_gf, only : green_stk, green_dip
   use rise_time, only : source, fourier_asym_cosine, realtr, fft
   use modelling_inputs, only : get_annealing_param
   implicit none
   integer, parameter :: nnsta_tele = 80


contains


   subroutine write_forward(slip, rake, rupt_time, tl, tr, strong, cgnss, body, surf, dart)
   !!
   !!  We save the forward solution given a kinematic model, for all specified 
   !!  data types.
   !!  
   implicit none
   integer ll_in, ll_out
   logical :: strong, cgnss, body, surf, dart
   real slip(nnxy, max_seg), rake(nnxy, max_seg), rupt_time(nnxy, max_seg), tr(nnxy, max_seg), &
   &  tl(nnxy, max_seg), erm, ermin
   complex z0
!
   z0 = cmplx(0.0, 0.0)
   erm = 0.0
   ll_in = 0
   ll_out = 0
   ermin = 1.0e+10
   write(*,*) dxs, dys
   if (strong) then
      call write_strong_motion_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
      ll_in = ll_out
   end if
   if (cgnss) then
      call write_cgnss_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
      ll_in = ll_out
   end if
   if (body) then
      call write_body_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
      ll_in = ll_out
   end if
   if (surf) then
      call write_surface_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
      ll_in = ll_out
   end if
   if (dart) then
      call write_dart_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
      ll_in = ll_out
   end if
!   select case (io_data)
!      case (1)
!         call write_cgnss_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!      case (2)
!         call write_strong_motion_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!      case (3)
!         call write_body_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!      case (4)
!         call write_strong_motion_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         lL_in = ll_out
!         call write_cgnss_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!      case (5)
!         call write_body_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_cgnss_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!      case (6)
!         call write_body_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_strong_motion_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!      case (7)
!         call write_strong_motion_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_surface_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!      case (8)
!         call write_cgnss_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_surface_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!      case (9)
!         call write_body_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_surface_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!      case (10)
!         call write_body_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_strong_motion_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = lL_out
!         call write_cgnss_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!      case (11)
!         call write_body_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_surface_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_cgnss_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!      case (12)
!         call write_body_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_surface_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_strong_motion_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!      case (13)
!         call write_strong_motion_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_surface_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         lL_in = ll_out
!         call write_cgnss_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!      case (14)
!         call write_body_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_strong_motion_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_cgnss_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!      case (15)
!         call write_body_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_strong_motion_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_surface_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!         ll_in = ll_out
!         call write_cgnss_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
!   endselect
   end subroutine write_forward
   
   
   subroutine write_strong_motion_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, isl, isr, ll, ir_max, int1, int2, int3, &
   &  jf, i, k, ll_s, i_s, ir, n_chan, ixs, iys, n_chan3
   real slip(nnxy, max_seg), rake(nnxy, max_seg), rupt_time(nnxy, max_seg), &
   &  tr(nnxy, max_seg), tl(nnxy, max_seg), &
   &  cr(inptd), cz(inptd), r, t1, t2, time, a, b, ww, float1, float2, float3,&
   &  dt, rake2, umax, df
   complex forward(inptd), z0, z
   complex :: source2(npth, mmsou, mmsou)
   character(len=6) sta_name(nnsta)
   character(len=3) component(nnsta), comp

   z0 = cmplx(0.0, 0.0)
   n_chan3 = 0
   open(9,file='Readlp.inf',status='old')
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*) ir_max, n_chan
   read(9,*)
   do ir = 1, ir_max
      read(9,*) int1, sta_name(ir), float1, float2, int2, component(ir), float3, int3
   end do
   close(9)
!
!       make the rise time function
!   
   dt = dt_channel(ll_in + 1)
   jf = 2**(lnpt-1)+1
   df = 1.0/(2**lnpt)/dt
   open(120,file='Source_function')
  
   do isl = 1, msou
      do isr = 1, msou
         cr(:) = 0.0
         cz(:) = 0.0
         t1 = ta0+(isl-1)*dta
         t2 = ta0+(isr-1)*dta
         r = t1+t2 
         k = int(r/dt+0.5)+1
         do i = 1, k
            time = (i-1)*dt
            if (time .lt. t1) then
               cr(i) = (1.0-cos(2*pi*time/(2*t1)))/r
            else
               cr(i) = (1.0+cos(2*pi*(time-t1)/(2*t2)))/r
            end if
         end do
         write(120,*)'k: ', k
         write(120,*)'(t1, t2): ', t1, t2
         write(120,*)'STF: ', (cr(i),i=1,512)
         CALL FFT(CR, CZ, LNPT,-1.)
         do i = 1, jf
            write(120,*)'freq: ', (i - 1) * df
            write(120,*)'F[STF]: ', (cmplx(cr(i),cz(i))) * dt
            call fourier_asym_cosine((i-1)*df, t1, t2, source2(i, isl, isr))
!                SOURCE(i, ISL, isr) = CMPLX(CR(i), CZ(i))*dt
            write(120,*)'F[STF]: ', source2(i,isl,isr)
         end do
      end do
   end do

   close(120)
!
!  end of rise time 
!       
   open(18,file='synm.str')
   k = 0
!       
!  set up the green function for every subfault
!  and calculate the initial value of objective function
!
   do ir = 1, n_chan
      comp = component(ir)
   
      ll_g = ir+ll_in
      do i = 1, npth
         forward(i) = z0
      end do
      ll = 0
      do i_s = 1, n_seg
         do iys = 1, nys_sub(i_s)
            do ixs = 1, nxs_sub(i_s)
               ll = ll+1  
               ll_s = (iys-1)*nxs_sub(i_s)+ixs             
               isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
               isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
               rake2 = rake(ll_s, i_s)*dpi
               a = sin(rake2)*slip(ll_s, i_s)
               b = cos(rake2)*slip(ll_s, i_s)
               do i = 1, nlen 
                  ww = -(i-1)*twopi*df*rupt_time(ll_s, i_s)
                  z = cmplx(cos(ww), sin(ww))
                  forward(i) = forward(i) &
               &  +(a*green_dip(i, ll_g, ll)+b*green_stk(i, ll_g, ll))*source2(i, isl, isr)*z
               end do
            end do
         end do
      end do

      do i = 1, jf
         if (i .le. jfmax) then
            cr(i) = real(forward(i))
            cz(i) = aimag(forward(i))
         else
            cr(i) = 0.0
            cz(i) = 0.0
         end if
      end do
       
      call realtr(cr, cz, lnpt)
      call fft(cr, cz, lnpt, 1.)
   
      umax = 0.0
    
      if (comp .eq.'HNZ') write(18,*)nlen,dt,sta_name(ir),'HNZ'
      if (comp .eq.'HNN') write(18,*)nlen,dt,sta_name(ir),'HNN'
      if (comp .eq.'HNE') write(18,*)nlen,dt,sta_name(ir),'HNE'
      do k = 1, nlen
         write(18,*) cr(k), cz(k)
      end do
   
   end do
   close(18)
   ll_out = ll_out+n_chan
   end subroutine write_strong_motion_forward
   
   
   subroutine write_cgnss_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, isl, isr, ll, n_chan3, int1, int2, int3, &
   &  jf, i, k, ll_s, i_s, ir, n_chan, ixs, iys, ir_max
   real slip(nnxy, max_seg), rake(nnxy, max_seg), rupt_time(nnxy, max_seg), &
   &  tr(nnxy, max_seg), tl(nnxy, max_seg), cr(inptd), cz(inptd), t1, t2, a, &
   &  b, ww, dt, rake2, df, float1, float2, float3
   complex forward(inptd), z0, z
   complex :: source2(npth, mmsou, mmsou)
   character(len=6) sta_name(nnsta)
   character(len=3) component(nnsta), comp

   z0 = cmplx(0.0, 0.0)
   n_chan3 = 0
   
   open(9,file='Readlp.cgnss',status='old')
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*) ir_max, n_chan
   read(9,*)
   do ir = 1, ir_max
      read(9,*) int1, sta_name(ir), float1, float2, int2, component(ir), float3, int3
   end do
   close(9)
!
!       make the rise time function
!     
   dt = dt_channel(ll_in + 1)
   jf = 2**(lnpt-2)+1
   df = 1.0/(2**lnpt)/dt
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
   open(18,file='synm.cgnss')
   k = 0
!       
!  set up the green function for every subfault
!  and calculate the initial value of objective function
!
   do ir = 1, n_chan
      comp = component(ir)
      ll_g = ir+ll_in
   
      do i = 1, npth
         cr(i) = 0.0
         cz(i) = 0.0
         forward(i) = z0
      end do
      ll = 0
      do i_s = 1, n_seg
         do iys = 1, nys_sub(i_s)
            do ixs = 1, nxs_sub(i_s)
               ll = ll+1  
               ll_s = (iys-1)*nxs_sub(i_s)+ixs             
               isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
               isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
               rake2 = rake(ll_s, i_s)*dpi
               a = sin(rake2)*slip(ll_s, i_s)
               b = cos(rake2)*slip(ll_s, i_s)
               do i = 1, nlen 
                  ww = -(i-1)*twopi*df*rupt_time(ll_s, i_s)
                  z = cmplx(cos(ww), sin(ww))
                  forward(i) = forward(i) &
               &  +(a*green_dip(i, ll_g, ll)+b*green_stk(i, ll_g, ll))*source2(i, isl, isr)*z
               end do
            end do
         end do
      end do
   
      do i = 1, jf
         cr(i) = real(forward(i))
         cz(i) = aimag(forward(i))
      end do
  
      call realtr(cr, cz, lnpt)
      call fft(cr, cz, lnpt, 1.)
   
      if (comp .eq.'LXZ') write(18,*)nlen,dt,sta_name(ir),'LXZ'
      if (comp .eq.'LXN') write(18,*)nlen,dt,sta_name(ir),'LXN'
      if (comp .eq.'LXE') write(18,*)nlen,dt,sta_name(ir),'LXE'
      do k = 1, nlen
         write(18,*) cr(k), cz(k)
      end do
   end do   
   close(18)
   ll_out = ll_out+n_chan
   end subroutine write_cgnss_forward


   subroutine write_body_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
   implicit none
   integer nnn, i_seg, nxy, nstaon, ir, ll_g, k, &
   &  ll, i_s, iys, jf, i, npxy, ll_s, kxy, ixs, &
   &  isl, isr, nl, llove(nnsta_tele), int1, int2, &
   &  nxys(max_seg), ll_in, ll_out, n_chan3
   real slip(nnxy, max_seg), rake(nnxy, max_seg), rupt_time(nnxy, max_seg), &
   &  tr(nnxy, max_seg), tl(nnxy, max_seg), t1, t2, &
   &  dt, df, ddelt, azim, w, cr(inptd), cz(inptd), sinal, cosal, &
   &  float1, float2, float3, float4, float5, float6, float7
   complex ::  z, z0, forward(npth)
   complex :: source2(npth, mmsou, mmsou)
!
   character(len=6) STNAME(nnsta_tele), string1, string2, string3
!
   open(9,file='Readlp.das',status='old')
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*) nstaon
   do ir = 1, nstaon
      read(9,*) int1, string1, string2, stname(ir), string3, &
              & float1, float2, float3, float3, float4, float5, int2, float6, float7, llove(ir)
   end do
   close(9)

   z0 = cmplx(0.0, 0.0)
   n_chan3 = 0
   nnn = 0
   do i_seg = 1, n_seg
      nxy = nxs_sub(i_seg)*nys_sub(i_seg)
      nxys(i_seg) = nxy
      nnn = nnn+nxy
   end do

   dt = dt_channel(ll_in + 1)
   jf = 2**(lnpt-1)+1
   df = 1.0/(2**lnpt)/dt
   do isl = 1, msou
      do isr = 1, msou
         t1 = ta0+(isl-1)*dta
         t2 = ta0+(isr-1)*dta
         if (t1 .lt. dt) t1 = dt
         if (t2 .lt. dt) t2 = dt
         do i = 1, 2*jfmax
            call fourier_asym_cosine((i-1)*df, t1, t2, source2(i, isl, isr))
         end do
      end do
   end do
!
!  End of Rise Time 
!
   OPEN(18,FILE = 'synm.tele')

   k = 0
!       
!  Now, we compute the synthetic seismographs
!
   npxy = nx_p*ny_p
   do ir = 1, nstaon
      ll_g = ll_in+ir

      do i = 1, npth
         forward(i) = z0
      end do
      LL = 0
      ddelt = 0.0
      do i_s = 1, n_seg
         kxy = 0
         do iys = 1, nys_sub(i_s)
            do ixs = 1, NXS_sub(i_s)
               kxy = kxy+1
               LL = LL+1
 
               azim = rake(kxy, i_s)*dpi
               sinal = sin(azim)*slip(kxy, i_s)
               cosal = cos(azim)*slip(kxy, i_s)
               ll_s = (iys-1)*nxs_sub(i_s)+ixs             
               isl = int((tl(kxy, i_s)-ta0)/dta+0.5)+1
               isr = int((tr(kxy, i_s)-ta0)/dta+0.5)+1
               do i = 1, 2*jfmax      
                  w = -(i-1)*twopi*df*rupt_time(kxy, i_s)
                  z = cmplx(cos(w), sin(w))
                  forward(i) = forward(i)&
               & +(sinal*green_dip(i, ll_g, ll)+cosal*green_stk(i, ll_g, ll))*source2(i, Isl, isr)*z
               end do
            end do
         end do
      end do
      do i = 1, jf
         if (i .le. jfmax) then
            cr(i) = real(forward(i))
            cz(i) = aimag(forward(i))
         else
            cr(i) = 0.0
            cz(i) = 0.0
         end if
      end do
      call realtr(cr, cz, lnpt)
      call fft(cr, cz, lnpt, 1.0)
      nl = 2**lnpt
      if (llove(ir) .eq. 0) then
         write(18,*)nl,dt,stname(ir),'P'
      else
         write(18,*)nl,dt,stname(ir),'SH'
      end if
      do i = 1, nl
         write(18,*) cr(i), cz(i)
      end do
   end do
!
   close(18)
   ll_out = ll_in+nstaon
   end subroutine write_body_waves_forward


   subroutine write_surface_waves_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ir_max, &
   &  ll_g, isl, isr, ll, jf, i, k, ll_s, i_s, ir, n_chan, &
   &  iys, ixs, io_up(nnsta), int1, int2

   real slip(nnxy, max_seg), rake(nnxy, max_seg), rupt_time(nnxy, max_seg), &
   &  tr(nnxy, max_seg), tl(nnxy, max_seg), &
   &  cr(inptd), cz(inptd), t1, t2, a, b, ww, dt, rake2, umax, df, float1, float2

   complex z0, forward(inptd), z
   complex :: source2(npth, mmsou, mmsou)

   character(len=6) sta_name(nnsta)

   z0 = cmplx(0.0, 0.0)

   open(9,file='Readlp.inf_low',status='old')
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*) ir_max, n_chan
   write(*,*) n_chan
   read(9,*)
   do ir = 1, ir_max
      read(9,*) int1, sta_name(ir), float1, float2, int2, io_up(ir)
   end do
   close(9)
!
! suppose the ni = u3e+11, then then moment of 1cm*1km^2 
! subfault is 3e+21. The gf is for Mo = 1e+20 dyne.cm
! the unit of surface wave is mm.
!
!       make the rise time function
!     
   dt = dt_channel(ll_in + 1)
   jf = 2**(lnpt-1)+1
   df = 1.0/(2**lnpt)/4.0
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
   open(18,file='synm.str_low')
   k = 0
!       
!  set up the green function for every subfault
!  and calculate the initial value of objective function
!
   do ir = 1, n_chan
      ll_g = ir+ll_in
   
      do i = 1, npth
         cr(i) = 0.0
         cz(i) = 0.0
         forward(i) = z0
      end do
      ll = 0
      do i_s = 1, n_seg
         do iys = 1, nys_sub(i_s)
            do ixs = 1, nxs_sub(i_s)
               ll = ll+1  
               ll_s = (iys-1)*nxs_sub(i_s)+ixs             
               isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
               isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
               rake2 = rake(ll_s, i_s)*dpi
               a = sin(rake2)*slip(ll_s, i_s)
               b = cos(rake2)*slip(ll_s, i_s)
               do i = 1, jfmax
                  ww = -(i-1)*twopi*df*rupt_time(ll_s, i_s)
                  z = cmplx(cos(ww), sin(ww))
                  forward(i) = forward(i) &
               & +(a*green_dip(i, ll_g, ll)+b*green_stk(i, ll_g, ll))*source2(i, isl, isr)*z
               end do
            end do
         end do
      end do
      do i = 1, jfmax
         cr(i) = real(forward(i))
         cz(i) = aimag(forward(i))
      end do
     
      call realtr(cr, cz, lnpt)
      call fft(cr, cz, lnpt, 1.0)
   
      umax = 0.0
   
      if (io_up(ir) .eq. 1) then
         write(18,*)nlen,dt,sta_name(ir),'P'
      else
         write(18,*)nlen,dt,sta_name(ir),'SH'
      end if
      do k = 1, nlen
         write(18,*) cr(k), cz(k)
      end do
   end do   
   close(18)
   ll_out = ll_out+n_chan

   end subroutine write_surface_waves_forward


   subroutine write_dart_forward(slip, rake, rupt_time, tl, tr, ll_in, ll_out)
   implicit none
   integer ll_in, ll_out, ll_g, isl, isr, ll, n_chan3, int1, int2, int3, &
   &  jf, i, k, ll_s, i_s, ir, n_chan, ixs, iys, ir_max
   real slip(nnxy, max_seg), rake(nnxy, max_seg), rupt_time(nnxy, max_seg), &
   &  tr(nnxy, max_seg), tl(nnxy, max_seg), cr(inptd), cz(inptd), t1, t2, a, &
   &  b, ww, dt, rake2, df, float1, float2, float3
   complex forward(inptd), z0, z
   complex :: source2(npth, mmsou, mmsou)
   character(len=6) sta_name(nnsta)

   z0 = cmplx(0.0, 0.0)
   n_chan3 = 0
   
   open(9,file='Readlp.dart',status='old')
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*)
   read(9,*) ir_max, n_chan
   read(9,*)
   do ir = 1, ir_max
      read(9,*) int1, sta_name(ir), float1, float2, int2, float3, int3
   end do
   close(9)
!
!       make the rise time function
!     
   dt = dt_channel(ll_in + 1)
   jf = 2**(lnpt-2)+1
   df = 1.0/(2**lnpt)/dt
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
   open(18,file='synm.dart')
   k = 0
!       
!  set up the green function for every subfault
!  and calculate the initial value of objective function
!
   do ir = 1, n_chan
      ll_g = ir+ll_in
   
      do i = 1, npth
         cr(i) = 0.0
         cz(i) = 0.0
         forward(i) = z0
      end do
      ll = 0
      do i_s = 1, n_seg
         do iys = 1, nys_sub(i_s)
            do ixs = 1, nxs_sub(i_s)
               ll = ll+1  
               ll_s = (iys-1)*nxs_sub(i_s)+ixs             
               isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
               isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
               rake2 = rake(ll_s, i_s)*dpi
               a = sin(rake2)*slip(ll_s, i_s)
               b = cos(rake2)*slip(ll_s, i_s)
               do i = 1, nlen 
                  ww = -(i-1)*twopi*df*rupt_time(ll_s, i_s)
                  z = cmplx(cos(ww), sin(ww))
                  forward(i) = forward(i) &
               &  +(a*green_dip(i, ll_g, ll)+b*green_stk(i, ll_g, ll))*source2(i, isl, isr)*z
               end do
            end do
         end do
      end do
   
      do i = 1, jf
         cr(i) = real(forward(i))
         cz(i) = aimag(forward(i))
      end do
  
      call realtr(cr, cz, lnpt)
      call fft(cr, cz, lnpt, 1.)
   
      write(18,*)nlen,dt,sta_name(ir), 'dart'
      do k = 1, nlen
         write(18,*) cr(k), cz(k)
      end do
   end do   
   close(18)
   ll_out = ll_out+n_chan
   end subroutine write_dart_forward


end module save_forward
        
