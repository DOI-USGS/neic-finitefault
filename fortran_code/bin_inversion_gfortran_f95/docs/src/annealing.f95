module annealing


   use constants, only : pi, nnsta, dpi, twopi, nnxy, max_seg, inptd, nt1, npth, nnxy_m
   use retrieve_gf, only : green_dip, green_stk
   use wavelets, only : wavelet_syn
   use wavelet_param, only : jfmax, lnpt 
   use rise_time, only : source
   use get_stations_data, only : dt_channel
   use random_gen, only : ran1, cauchy
   use misfit_eval, only : misfit
   use modelling_inputs, only : smooth_moment, smooth_slip, smooth_time, io_re, cm_point
   use model_parameters, only : n_seg, ta0, dta, msou, dxs, dys, nxs0, nys0, nx_p, ny_p, v_min, v_max, &
                            &   tbl, tbr, v_ref, nxs_sub, nys_sub, time_max, time_min, cniu, &
                            &   slip0, rake0, rupt_time0, tl0, tr0, c_depth, beg, dp, np, dip_seg, &
                            &   stk_seg, delay_seg
   use regularization, only : lap, tlap, define_slip_field, modify_slip_field
   use static_data, only : static_synthetic, static_remove_subfault, &
                       &   static_modify_subfault, static_add_subfault
   implicit none
   real :: coef_moment, coef_slip, coef_st, coef_time
   real :: emin, ermin, min_dt, area
   integer :: n_chan, nxys(max_seg), nnn
   integer, parameter :: double = kind(1.d0)
   integer, private :: threads
   integer, parameter, private :: max_move=50, accept_max=5


contains


   subroutine n_threads(auto)
   use omp_lib
   implicit none
   logical :: auto
   threads = 4
   if (auto) threads = 3
   call omp_set_num_threads(threads)
   end subroutine n_threads
   
   
   subroutine initial_model(slip, rake, rupt_time, tl, tr)
   !!
   !!  We define initial model for the annealing method.
   !!
   implicit none
   real :: slip(nnxy, max_seg), rake(nnxy, max_seg), rupt_time(nnxy, max_seg), &
         & tl(nnxy, max_seg), tr(nnxy, max_seg)
   real :: x(nt1)
   integer :: nxy, i, i_g, i_seg, io_v_d, k, kp, npa
   io_v_d = 0
   npa = 0
   do i_seg = 1, n_seg
      npa = npa + 4*nxs_sub(i_seg)*nys_sub(i_seg)
   end do
   do k = 1, npa
      x(k) = beg(k)+int(ran1()*np(k)-1.0)*dp(k)
   end do
   if (io_re .eq. 0) then
      k = 0
      do i_seg = 1, n_seg
         nxy = nxs_sub(i_seg)*nys_sub(i_seg)
         do i = 1, nxy
            k = k + 1
            slip(i, i_seg) = x(k)
            k = k + 1
            rake(i, i_seg) = x(k)
            k = k + 2
            tl(i, i_seg) = (ta0+int(ran1()*msou)*dta)
            tr(i, i_seg) = (ta0+int(ran1()*msou)*dta)
            rupt_time(i, i_seg) = ran1()* &
     &  (time_max(i, i_seg)-time_min(i, i_seg))+time_min(i, i_seg)
         end do
      end do 
      open(12,file='Fault.ini')
      write(12,*) nxs0, nys0, c_depth
      write(12,*) n_seg, dxs, dys, nx_p, ny_p, v_min, v_max, tbl, tbr
      write(12,*) ta0, dta, msou, v_ref, io_v_d
      do i_g = 1, n_seg
         nxys(i_g) = nxs_sub(i_g)*nys_sub(i_g)
         write(12,*) i_g, dip_seg(i_g), stk_seg(i_g)
         write(12,*) nxs_sub(i_g), nys_sub(i_g), delay_seg(i_g)
         do kp = 1, nxys(i_g)
            write(12,*) slip(kp, i_g), rake(kp, i_g), rupt_time(kp, i_g), tl(kp, i_g), tr(kp, i_g)
         end do
      end do
      close(12)
   else
      slip(:, :) = slip0(:, :)
      rake(:, :) = rake0(:, :)
      rupt_time(:, :) = rupt_time0(:, :)
      tl(:, :) = tl0(:, :)
      tr(:, :) = tr0(:, :)
   end if
   end subroutine initial_model


   subroutine initial_regularization(slip, rake, rupt_time, tl, tr, static)
   !!
   !!  We compute regularization coefficients, and store them in memory.
   !!
   implicit none
   real :: slip(nnxy, max_seg), rake(nnxy, max_seg), rupt_time(nnxy, max_seg)
   real :: tl(nnxy, max_seg), tr(nnxy, max_seg)
   real amp, cmoment, derr, dt, er, er0, err_s, err_static, &
      & err_time, forward_real(inptd, nnsta), a, b, &
      & forward_imag(inptd, nnsta), cr(inptd), cz(inptd), forward2(inptd)
   real :: ex, rake2, df, df_ref, dsum, kahan_y, kahan_t, kahan_c
   real*8 :: omega, ersum2
   integer :: i, i_s, ir, isl, isr, ixs, iys, jf, k, ll, ll_s 
   complex*16 :: z0, forward(npth), z, z1
   logical :: static

   z0 = cmplx(0.d0, 0.d0, double)
   min_dt = 10
   do i = 1, nnsta
      if (dt_channel(i) .gt. 1.e-4) min_dt = min(min_dt, dt_channel(i))
   end do
   n_chan = 0
   do i = 1, nnsta
      if (dt_channel(i) .gt. 1.e-4) n_chan = n_chan + 1
   end do
   jf = 2**(lnpt-1)+1
!
! Compute synthetics given current fault model
!
   area = dxs*dys*(1.e+10)
   df_ref = 1.0/(2.0**lnpt)
   ersum2 = 0.0
   do ir = 1, n_chan
      df = df_ref/dt_channel(ir)
      dt = dt_channel(ir)
      do i = 1, npth
         forward(i) = z0
      end do
      ll = 0
      do i_s = 1, n_seg
         do iys = 1, nys_sub(i_s)
            do ixs = 1, nxs_sub(i_s)
               ll = ll+1
               ll_s = (iys-1)*nxs_sub(i_s)+ixs
               rake2 = rake(ll_s, i_s)*dpi
               a = sin(rake2)*slip(ll_s, i_s)
               b = cos(rake2)*slip(ll_s, i_s)
               isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
               isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
               omega = -twopi*df*rupt_time(ll_s, i_s)
               z1 = cmplx(cos(omega), sin(omega), double)
               z = cmplx(1.d0, 0.d0, double)
               do i = 1, jfmax
                  forward(i) = forward(i) &
               &  +(a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll))*source(i, ir, isl, isr)*z
                  z = z*z1
               end do
            end do
         end do
      end do

      do i = 1, npth
         forward_real(i, ir) = real(forward(i))
         forward_imag(i, ir) = aimag(forward(i))
         cr(i) = forward_real(i, ir)
         cz(i) = forward_imag(i, ir)
      end do
      call wavelet_syn(cr, cz, forward2)
      call misfit(ir, forward2, ex)
      ersum2 = ersum2 + ex
   end do

   amp = 1.0
   dsum = 0.0
   kahan_y = 0.0
   kahan_t = 0.0
   kahan_c = 0.0
   do i_s = 1, n_seg
      do k = 1, nxys(i_s)
         kahan_y = slip(k, i_s)*cniu(k, i_s)-kahan_c 
         kahan_t = dsum+kahan_y
         kahan_c = (kahan_t-dsum)-kahan_y
         dsum = kahan_t
!         dsum = dsum+slip(k, i_s)*cniu(k, i_s)    ! we may need to increase numerical precision
      end do
   end do
   
   cmoment = dsum*area
   derr = (cmoment/cm_point)
   er0 = ersum2
   coef_moment = smooth_moment
   coef_moment = min(coef_moment, 1.0)
   call define_slip_field(slip, rake)
   call lap(err_s)
   coef_slip = smooth_slip*er0/(err_s*amp)
   call tlap(rupt_time, err_time)
   if (err_time .lt. 1.0e-9) err_time = 1.0
   coef_time = smooth_time*er0/(err_time*amp)
 
   coef_st = 0.0
   if (static) then
      call static_synthetic(slip, rake, nxys, err_static)
      coef_st = er0/(err_static*amp)
   end if
   coef_slip = min(0.003, coef_slip)
   er = ersum2+coef_moment*derr+coef_slip*err_s*amp
   er = er+coef_time*err_time+coef_st*err_static
   nnn = 0
   do i_s = 1, n_seg
      nxys(i_s) = nys_sub(i_s)*nxs_sub(i_s)
      nnn = nnn+nxys(i_s)
   end do
   write(*,*)''
   write(*,*)'averaged misfit error', er0
   write(*,*)'moment error', derr
   write(*,*)'slip smoothness penalization', err_s
   write(*,*)'time smoothness penalization', err_time
   if (static) write(*,*)'static data penalization', err_static
   write(*,*)'total moment of the inversion', cmoment
   write(*,*)''
   write(*,*)'moment error coefficient', coef_moment
   write(*,*)'slip smoothness penalization coefficient', coef_slip
   write(*,*)'time smoothness penalization coefficient', coef_time
   if (static) write(*,*)'static data penalization coefficient', coef_st
   write(*,*)''
   write(*,*)'Amount of variables: ', 5 * nnn
   write(*,*)'Amount of data values: ', n_chan * jfmax
   emin = er
   ermin = er
   end subroutine initial_regularization


   subroutine annealing_iter(slip, rake, rupt_time, tl, tr, er, t)
   !!
   !!  We compute an iteration of the annealing method.
   !!  Non-parallelized version.
   !!
   implicit none
   integer isl, isr, nn_sub, nsub(nnxy_m), n_accept, &
   & nbb, i, k, npb, nn, nran, ll_s, i_s, i_ss, ir, ll, iys, &
   & ixs, i_move, n_total
   real slip(nnxy, max_seg), rake(nnxy, max_seg), rupt_time(nnxy, max_seg), &
   & tr(nnxy, max_seg), tl(nnxy, max_seg), er, t, &
   & forward_real(npth, nnsta), forward_imag(npth, nnsta), duse, ause, vuse, &
   & de, rand, c, aux, dpb, amp, derr, erm, &
   & cmoment, d_sub, a_sub, err_s, a, b, dsum, &
   & err_time, t_save, a_save, d_save, x, kahan_y, kahan_t, kahan_c, &
   & l_save, r_save, cr(inptd), cz(inptd), forward2(inptd), &
   & slip_beg, slip_max, slip_end, angle_beg, angle_end, angle_max, &
   & rupt_beg, rupt_end, rupt_max, rise_time_beg, rise_time_end, rise_time_max
   real :: df, df_ref, rake2, ex
   real*8 :: ersum2, omega
   complex :: green_subf
   complex*16 :: z, z1, z0, forward(npth)
!
   z0 = cmplx(0.d0, 0.d0, double)
   erm = 0.0
!
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Here, we compute the value of the objective function, 
!  using the input kinematic model.
!
   ll = 0
   df_ref = 1.0/(2.0**lnpt)

   do ir = 1, n_chan
      df = df_ref/dt_channel(ir)
      forward(:) = z0
      ll = 0
      do i_s = 1, n_seg
         do iys = 1, nys_sub(i_s)
            do ixs = 1, nxs_sub(i_s)
               ll = ll+1
               ll_s = (iys-1)*nxs_sub(i_s)+ixs
               rake2 = rake(ll_s, i_s)*dpi
               a = sin(rake2)*slip(ll_s, i_s)
               b = cos(rake2)*slip(ll_s, i_s)
               isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
               isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
               omega = -twopi*df*rupt_time(ll_s, i_s)
               z1 = cmplx(cos(omega), sin(omega), double)
               z = cmplx(1.d0, 0.d0, double)
               do i = 1, jfmax
                  forward(i) = forward(i) &
                  & +(a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll))*source(i, ir, isl, isr)*z
                  z = z*z1
               end do
            end do
         end do
      end do

      do i = 1, npth
         forward_real(i, ir) = real(forward(i))
         forward_imag(i, ir) = aimag(forward(i))
      end do
   end do

   dsum = 0.0
   kahan_y = 0.0
   kahan_t = 0.0
   kahan_c = 0.0
   do i_s = 1, n_seg
      do k = 1, nxys(i_s)
         kahan_y = slip(k, i_s)*cniu(k, i_s)-kahan_c 
         kahan_t = dsum+kahan_y
         kahan_c = (kahan_t-dsum)-kahan_y
         dsum = kahan_t
!         dsum = dsum+slip(k, i_s)*cniu(k, i_s)    ! we may need to increase numerical precision
      end do
   end do

   call define_slip_field(slip, rake)
   call tlap(rupt_time, err_time)
!
!  An iteration of the simulated annealing algorithm, over each subfault.
!
   ll = 0
!       begin to perturb       
!
   do k = 1, nnn
      nsub(k) = k
   end do

   do k = 1, nnn-1
      nran = k
      do while (nran .eq. k .or. nran .gt. nnn)
         x = ran1()
         nran = int(x*(nnn-k)+k+1)
      end do
      nbb = nsub(nran)
      nn = nsub(k)
      nsub(k) = nbb
      nsub(nran) = nn
   end do

   do k = 1, nnn
      ll = nsub(k)
      if (ll .gt. nnn) stop
      n_total = 0
      do i_ss = 1, n_seg
         n_total = nxys(i_ss)+n_total
         if (ll .le. n_total) then
            i_s = i_ss
            ll_s = ll
            exit
         end if
      end do
      do i_ss = 1, i_s-1
         ll_s = ll_s-nxys(i_ss)
      end do
      isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
      isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
      rake2 = rake(ll_s, i_s)*dpi
      a = sin(rake2)*slip(ll_s, i_s)
      b = cos(rake2)*slip(ll_s, i_s)
      a_sub = rake(ll_s, i_s)
!
!  make up unchange graph
!
      do ir = 1, n_chan
         df = df_ref/dt_channel(ir)
         omega = -twopi*df*rupt_time(ll_s, i_s)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, jfmax
            green_subf = (a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll))* &
            &    source(i, ir, isl, isr)*z
            forward_real(i, ir) = forward_real(i, ir)-real(green_subf)
            forward_imag(i, ir) = forward_imag(i, ir)-aimag(green_subf)
            z = z*z1    ! we may need to increase numerical precision
         end do
      end do
      dsum = dsum-slip(ll_s, i_s)*cniu(ll_s, i_s)
      nn_sub = ll
!  
      n_accept = 0
      npb = np(4*(ll-1)+1)
      if (npb .lt. 2) exit
!
!  slip extreme values
!
      npb = np(4*(ll-1)+1)
      dpb = dp(4*(ll-1)+1)
      slip_beg = beg(4*(ll-1)+1)
      slip_max = (npb-1)*dpb
      slip_end = slip_beg+slip_max
!
!  rake extreme values
!  
      npb = np(4*(ll-1)+2)
      dpb = dp(4*(ll-1)+2)
      angle_beg = beg(4*(ll-1)+2)
      angle_max = (npb-1)*dpb
      angle_end = angle_beg+angle_max
! 
!  rupture time extreme values.
!
      rupt_beg = time_min(ll_s, i_s)
      rupt_end = time_max(ll_s, i_s)
      rupt_max = rupt_end-rupt_beg
      write(*,*)rupt_beg, rupt_end, rupt_max
!
!  rise time parameters extreme values
!  
      rise_time_beg = ta0
      rise_time_end = ta0+(msou-1)*dta
      rise_time_max = (msou-1)*dta
      do i_move = 1, max_move
!
!       Save values before the perturbation
!
         t_save = rupt_time(ll_s, i_s)
         d_save = slip(ll_s, i_s)
         a_save = rake(ll_s, i_s)
         l_save = tl(ll_s, i_s)
         r_save = tr(ll_s, i_s)
!
!  Perturb the slip
!
         duse = slip_beg - 1.
         do while ((duse .le. slip_beg) .or. (duse .ge. slip_end))
            call cauchy(t, c)                           
            duse = slip(ll_s, i_s)+c*slip_max
         end do
!
!  Perturb the rake
!
         ause = angle_beg - 1.
         do while ((ause .lt. angle_beg) .or. (ause .gt. angle_end))
            call cauchy(t, c)                          
            ause = rake(ll_s, i_s)+c*angle_max
         end do
! 
!  Perturb the rupture time.
!
         vuse = rupt_beg - 1.
         if (rupt_max .gt. min_dt) then
            do while ((vuse .lt. rupt_beg) .or. (vuse .gt. rupt_end))
               call cauchy(t, c)                       
               vuse = min_dt*int((rupt_time(ll_s, i_s)+c*rupt_max)/min_dt+0.5)     
            end do
         else
            vuse = rupt_beg
         end if             
!
!  Perturb rise time parameters
!  
         isl = 0
         do while (isl .lt. 1 .or. isl .gt. msou)
            call cauchy(t, c)                         
            isl = int((tl(ll_s, i_s)+c*rise_time_max-ta0)/dta+0.5)+1
         end do   
         isr = 0
         do while (isr .lt. 1 .or. isr .gt. msou)
            call cauchy(t, c)                         
            isr = int((tr(ll_s, i_s)+c*rise_time_max-ta0)/dta+0.5)+1
         end do
         
         rake2 = ause*dpi
         a = duse*sin(rake2)
         b = duse*cos(rake2)
         ersum2 = 0.d0
         do ir = 1, n_chan
            df = df_ref/dt_channel(ir)
            omega = -twopi*df*vuse
            z1 = cmplx(cos(omega), sin(omega), double) 
            z = cmplx(1.d0, 0.d0, double)
            do i = 1, jfmax
!               omega = -twopi_0*df*(i-1)*vuse
!               z = cmplx(cos(omega), sin(omega)) 
               green_subf = (a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll))* &
               &    source(i, ir, isl, isr)*z
               cr(i) = real(green_subf)+forward_real(i, ir)
               cz(i) = aimag(green_subf)+forward_imag(i, ir)
               z = z*z1    ! we may need to increase numerical precision
            end do
            call wavelet_syn(cr, cz, forward2)
            call misfit(ir, forward2, ex)     
            ersum2 = ersum2 + ex    ! we may need to increase numerical precision
         end do
         dsum = dsum+duse*cniu(ll_s, i_s)
         cmoment = dsum*area
         derr = (cmoment/cm_point)
         amp = 1.0
         d_sub = duse
         a_sub = ause
         call modify_slip_field(nn_sub, d_sub, a_sub)
         call lap(err_s)
         rupt_time(ll_s, i_s) = vuse
         call tlap(rupt_time, err_time)

         erm = ersum2+derr*coef_moment+amp*err_s*coef_slip
         erm = erm+coef_time*err_time
         dsum = dsum-duse*cniu(ll_s, i_s)
         de = erm-emin
!  
!  Now, we update the kinematic model.
!  
         rand = ran1()
         aux = exp(-de/t)
         if (aux .gt. rand) then
            emin = erm
            slip(ll_s, i_s) = duse
            rake(ll_s, i_s) = ause
            rupt_time(ll_s, i_s) = vuse
            tl(ll_s, i_s) = (isl-1)*dta+ta0
            tr(ll_s, i_s) = (isr-1)*dta+ta0
            n_accept = n_accept+1
         else
            rupt_time(ll_s, i_s) = t_save
            tl(ll_s, i_s) = l_save
            tr(ll_s, i_s) = r_save
            slip(ll_s, i_s) = d_save
            rake(ll_s, i_s) = a_save
         end if
         ermin = min(ermin, erm)
         if (n_accept .gt. accept_max) exit
      end do
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  finish the perturbation of subevent (i_s, ll_s)
!
      rake2 = rake(ll_s, i_s)*dpi
      a = sin(rake2)*slip(ll_s, i_s)
      b = cos(rake2)*slip(ll_s, i_s)
      isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
      isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
      dsum = dsum+slip(ll_s, i_s)*cniu(ll_s, i_s)
      d_sub = slip(ll_s, i_s)
      a_sub = rake(ll_s, i_s)
      nn_sub = ll
      call modify_slip_field(nn_sub, d_sub, a_sub)
      call tlap(rupt_time, err_time)
      do ir = 1, n_chan
         df = df_ref/dt_channel(ir)
         omega = -twopi*df*rupt_time(ll_s, i_s)
         z1 = cmplx(cos(omega), sin(omega), double)   
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, jfmax
            green_subf = (a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll)) &
            &   *source(i, ir, isl, isr)*z
            forward_real(i, ir) = forward_real(i, ir)+real(green_subf)
            forward_imag(i, ir) = forward_imag(i, ir)+aimag(green_subf)
            z = z*z1
         end do
      end do

   end do
   write(*,*) ermin

   end subroutine annealing_iter 
   
   
   subroutine annealing_iter2(slip, rake, rupt_time, tl, tr, er, t)
   !!
   !!  We compute an iteration of the annealing method, adding static GNSS data.
   !!  Non-parallelized version.
   !!
   implicit none
   integer isl, isr, nn_sub, nsub(nnxy_m), n_accept, &
   & nbb, i, k, npb, nn, nran, ll_s, i_s, i_ss, ir, ll, iys, &
   & ixs, i_move, n_total
   real slip(nnxy, max_seg), rake(nnxy, max_seg), rupt_time(nnxy, max_seg), &
   & tr(nnxy, max_seg), tl(nnxy, max_seg), er, t, &
   & forward_real(npth, nnsta), forward_imag(npth, nnsta), duse, ause, vuse, &
   & de, rand, c, aux, dpb, amp, derr, erm, err_static, &
   & cmoment, d_sub, a_sub, err_s, a, b, dsum, &
   & err_time, t_save, a_save, d_save, x, kahan_y, kahan_t, kahan_c, &
   & l_save, r_save, cr(inptd), cz(inptd), forward2(inptd), &
   & slip_beg, slip_max, slip_end, angle_beg, angle_end, angle_max, &
   & rupt_beg, rupt_end, rupt_max, rise_time_beg, rise_time_end, rise_time_max
   real :: df, df_ref, rake2, ex
   real*8 :: omega, ersum2
   complex :: green_subf
   complex*16 :: z, z1, forward(npth), z0
!
   z0 = cmplx(0.d0, 0.d0, double)
   erm = 0.0
!
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Here, we compute the value of the objective function, 
!  using the input kinematic model.
!
   ll = 0
   df_ref = 1.0/(2.0**lnpt)

   do ir = 1, n_chan
      df = df_ref/dt_channel(ir)
      forward(:) = z0
      ll = 0
      do i_s = 1, n_seg
         do iys = 1, nys_sub(i_s)
            do ixs = 1, nxs_sub(i_s)
               ll = ll+1
               ll_s = (iys-1)*nxs_sub(i_s)+ixs
               rake2 = rake(ll_s, i_s)*dpi
               a = sin(rake2)*slip(ll_s, i_s)
               b = cos(rake2)*slip(ll_s, i_s)
               isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
               isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
               omega = -twopi*df*rupt_time(ll_s, i_s)
               z1 = cmplx(cos(omega), sin(omega), double)
               z = cmplx(1.d0, 0.d0, double)
               do i = 1, jfmax
                  forward(i) = forward(i) &
                  & +(a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll))*source(i, ir, isl, isr)*z
                  z = z*z1    ! we may need to increase numerical precision
               end do
            end do
         end do
      end do

      do i = 1, npth
         forward_real(i, ir) = real(forward(i))
         forward_imag(i, ir) = aimag(forward(i))
      end do
   end do
   call static_synthetic(slip, rake, nxys, err_static)

   dsum = 0.0
   kahan_y = 0.0
   kahan_t = 0.0
   kahan_c = 0.0
   do i_s = 1, n_seg
      do k = 1, nxys(i_s)
         kahan_y = slip(k, i_s)*cniu(k, i_s)-kahan_c 
         kahan_t = dsum+kahan_y
         kahan_c = (kahan_t-dsum)-kahan_y
         dsum = kahan_t
!         dsum = dsum+slip(k, i_s)*cniu(k, i_s)    ! we may need to increase numerical precision
      end do
   end do

   call define_slip_field(slip, rake)
   call tlap(rupt_time, err_time)
!
!  An iteration of the simulated annealing algorithm, over each subfault.
!
   ll = 0
!       begin to perturb       
!
   do k = 1, nnn
      nsub(k) = k
   end do

   do k = 1, nnn-1
      nran = k
      do while (nran .eq. k .or. nran .gt. nnn)
         x = ran1()
         nran = int(x*(nnn-k)+k+1)
      end do
      nbb = nsub(nran)
      nn = nsub(k)
      nsub(k) = nbb
      nsub(nran) = nn
   end do

   do k = 1, nnn
      ll = nsub(k)
      if (ll .gt. nnn) stop
      n_total = 0
      do i_ss = 1, n_seg
         n_total = nxys(i_ss)+n_total
         if (ll .le. n_total) then
            i_s = i_ss
            ll_s = ll
            exit
         end if
      end do
      do i_ss = 1, i_s-1
         ll_s = ll_s-nxys(i_ss)
      end do
      isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
      isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
      rake2 = rake(ll_s, i_s)*dpi
      a = sin(rake2)*slip(ll_s, i_s)
      b = cos(rake2)*slip(ll_s, i_s)
      a_sub = rake(ll_s, i_s)
!
!  make up unchange graph
!
      do ir = 1, n_chan
         df = df_ref/dt_channel(ir)
         omega = -twopi*df*rupt_time(ll_s, i_s)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, jfmax
            green_subf = (a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll))* &
            &    source(i, ir, isl, isr)*z
            forward_real(i, ir) = forward_real(i, ir)-real(green_subf)
            forward_imag(i, ir) = forward_imag(i, ir)-aimag(green_subf)
            z = z*z1    ! we may need to increase numerical precision
         end do
      end do
      call static_remove_subfault(slip, rake, i_s, ll_s)
      dsum = dsum-slip(ll_s, i_s)*cniu(ll_s, i_s)
      nn_sub = ll
!  
      n_accept = 0
      npb = np(4*(ll-1)+1)
      if (npb .lt. 2) exit
!
!  slip extreme values
!
      npb = np(4*(ll-1)+1)
      dpb = dp(4*(ll-1)+1)
      slip_beg = beg(4*(ll-1)+1)
      slip_max = (npb-1)*dpb
      slip_end = slip_beg+slip_max
!
!  rake extreme values
!  
      npb = np(4*(ll-1)+2)
      dpb = dp(4*(ll-1)+2)
      angle_beg = beg(4*(ll-1)+2)
      angle_max = (npb-1)*dpb
      angle_end = angle_beg+angle_max
! 
!  rupture time extreme values.
!
      rupt_beg = time_min(ll_s, i_s)
      rupt_end = time_max(ll_s, i_s)
      rupt_max = rupt_end-rupt_beg
!
!  rise time parameters extreme values
!  
      rise_time_beg = ta0
      rise_time_end = ta0+(msou-1)*dta
      rise_time_max = (msou-1)*dta
      do i_move = 1, max_move
!
!       Save values before the perturbation
!
         t_save = rupt_time(ll_s, i_s)
         d_save = slip(ll_s, i_s)
         a_save = rake(ll_s, i_s)
         l_save = tl(ll_s, i_s)
         r_save = tr(ll_s, i_s)
!
!  Perturb the slip
!
         duse = slip_beg - 1.
         do while ((duse .le. slip_beg) .or. (duse .ge. slip_end))
            call cauchy(t, c)                           
            duse = d_save+c*slip_max
         end do
!
!  Perturb the rake
!
         ause = angle_beg - 1.
         do while ((ause .lt. angle_beg) .or. (ause .gt. angle_end))
            call cauchy(t, c)                          
            ause = a_save+c*angle_max
         end do
! 
!  Perturb the rupture time.
!
         vuse = rupt_beg - 1.
         if (rupt_max .gt. min_dt) then
            do while ((vuse .lt. rupt_beg) .or. (vuse .gt. rupt_end))
               call cauchy(t, c)                       
               vuse = min_dt*int((t_save+c*rupt_max)/min_dt+0.5)     
            end do
         else
            vuse = rupt_beg
         end if             
!
!  Perturb rise time parameters
!  
         isl = 0
         do while (isl .lt. 1 .or. isl .gt. msou)
            call cauchy(t, c)                         
            isl = int((l_save+c*rise_time_max-ta0)/dta+0.5)+1
         end do   
         isr = 0
         do while (isr .lt. 1 .or. isr .gt. msou)
            call cauchy(t, c)                         
            isr = int((r_save+c*rise_time_max-ta0)/dta+0.5)+1
         end do
         
         rake2 = ause*dpi
         a = duse*sin(rake2)
         b = duse*cos(rake2)
         ersum2 = 0.d0
         do ir = 1, n_chan
            df = df_ref/dt_channel(ir)
            omega = -twopi*df*vuse
            z1 = cmplx(cos(omega), sin(omega), double) 
            z = cmplx(1.d0, 0.d0, double)
            do i = 1, jfmax
!               omega = -twopi_0*df*(i-1)*vuse
!               z = cmplx(cos(omega), sin(omega)) 
               green_subf = (a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll))* &
               &    source(i, ir, isl, isr)*z
               cr(i) = real(green_subf)+forward_real(i, ir)
               cz(i) = aimag(green_subf)+forward_imag(i, ir)
               z = z*z1    ! we may need to increase numerical precision
            end do
            call wavelet_syn(cr, cz, forward2)
            call misfit(ir, forward2, ex)     
            ersum2 = ersum2 + ex    ! we may need to increase numerical precision
         end do
         call static_modify_subfault(duse, ause, i_s, ll_s, err_static)
         dsum = dsum+duse*cniu(ll_s, i_s)
         cmoment = dsum*area
         derr = (cmoment/cm_point)
         amp = 1.0
         d_sub = duse
         a_sub = ause
         call modify_slip_field(nn_sub, d_sub, a_sub)
         call lap(err_s)
         rupt_time(ll_s, i_s) = vuse
         call tlap(rupt_time, err_time)

         erm = ersum2+derr*coef_moment+amp*err_s*coef_slip
         erm = erm+coef_time*err_time+coef_st*err_static
         dsum = dsum-duse*cniu(ll_s, i_s)
         de = erm-emin
!  
!  Now, we update the kinematic model.
!  
         rand = ran1()
         aux = exp(-de/t)
         if (k .eq. 1) then 
            write(*,*) i_move, ersum2, err_s, err_static, err_time
            write(*,*) n_accept
         end if
         if (aux .gt. rand) then
            emin = erm
            slip(ll_s, i_s) = duse
            rake(ll_s, i_s) = ause
            rupt_time(ll_s, i_s) = vuse
            tl(ll_s, i_s) = (isl-1)*dta+ta0
            tr(ll_s, i_s) = (isr-1)*dta+ta0
            n_accept = n_accept+1
         else
            rupt_time(ll_s, i_s) = t_save
            tl(ll_s, i_s) = l_save
            tr(ll_s, i_s) = r_save
            slip(ll_s, i_s) = d_save
            rake(ll_s, i_s) = a_save
         end if
         ermin = min(ermin, erm)
         if (n_accept .gt. accept_max) exit
      end do
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  finish the perturbation of subevent (i_s, ll_s)
!
      rake2 = rake(ll_s, i_s)*dpi
      a = sin(rake2)*slip(ll_s, i_s)
      b = cos(rake2)*slip(ll_s, i_s)
      isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
      isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
      dsum = dsum+slip(ll_s, i_s)*cniu(ll_s, i_s)
      d_sub = slip(ll_s, i_s)
      a_sub = rake(ll_s, i_s)
      nn_sub = ll
      call modify_slip_field(nn_sub, d_sub, a_sub)
      call tlap(rupt_time, err_time)
      do ir = 1, n_chan
         df = df_ref/dt_channel(ir)
         omega = -twopi*df*rupt_time(ll_s, i_s)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, jfmax
!            omega = -twopi_0*df*(i-1)*rupt_time(ll_s, i_s)
!            z = cmplx(cos(omega), sin(omega))   
            green_subf = (a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll)) &
            &   *source(i, ir, isl, isr)*z
            forward_real(i, ir) = forward_real(i, ir)+real(green_subf)
            forward_imag(i, ir) = forward_imag(i, ir)+aimag(green_subf)
            z = z*z1
         end do
      end do
      call static_add_subfault(slip, rake, i_s, ll_s, err_static)

   end do
   write(*,*) ermin
   end subroutine annealing_iter2 
   
   
   subroutine annealing_iter3(slip, rake, rupt_time, tl, tr, er, t)
   !!
   !!  We compute an iteration of the annealing method.
   !!  Parallelized version.
   !!
   use omp_lib
   implicit none
   integer isl, isr, nn_sub, nsub(nnxy_m), n_accept, &
   & nbb, i, k, npb, nn, nran, ll_s, i_s, i_ss, ir, ll, iys, &
   & ixs, i_move, n_total
   real slip(nnxy, max_seg), rake(nnxy, max_seg), rupt_time(nnxy, max_seg), &
   & tr(nnxy, max_seg), tl(nnxy, max_seg), er, t, &
   & forward_real(npth, nnsta), forward_imag(npth, nnsta), duse, ause, vuse, &
   & de, rand, c, aux, dpb, amp, derr, erm, err_static, &
   & cmoment, d_sub, a_sub, err_s, a, b, dsum, &
   & err_time, t_save, a_save, d_save, x, kahan_y, kahan_t, kahan_c, &
   & l_save, r_save, cr(inptd), cz(inptd), forward2(inptd), &
   & slip_beg, slip_max, slip_end, angle_beg, angle_end, angle_max, &
   & rupt_beg, rupt_end, rupt_max, rise_time_beg, rise_time_end, rise_time_max
   real :: df_ref, df, rake2, ex
   real*8 :: omega, ersum2
   complex :: green_subf
   complex*16 :: z, z1, forward(npth), z0
!
   z0 = cmplx(0.d0, 0.d0, double)
   erm = 0.0
!
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Here, we compute the value of the objective function, 
!  using the input kinematic model.
!
   ll = 0
   df_ref = 1.0/(2.0**lnpt)

   do ir = 1, n_chan
      df = df_ref/dt_channel(ir)
      forward(:) = z0
      ll = 0
      do i_s = 1, n_seg
         do iys = 1, nys_sub(i_s)
            do ixs = 1, nxs_sub(i_s)
               ll = ll+1
               ll_s = (iys-1)*nxs_sub(i_s)+ixs
               rake2 = rake(ll_s, i_s)*dpi
               a = sin(rake2)*slip(ll_s, i_s)
               b = cos(rake2)*slip(ll_s, i_s)
               isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
               isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
               omega = -twopi*df*rupt_time(ll_s, i_s)
               z1 = cmplx(cos(omega), sin(omega), double)
               z = cmplx(1.d0, 0.d0, double)
               do i = 1, jfmax
                  forward(i) = forward(i) &
                  & +(a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll))*source(i, ir, isl, isr)*z
                  z = z*z1    ! we may need to increase numerical precision
               end do
            end do
         end do
      end do

      do i = 1, npth
         forward_real(i, ir) = real(forward(i))
         forward_imag(i, ir) = aimag(forward(i))
      end do
   end do

   dsum = 0.0
   kahan_y = 0.0
   kahan_t = 0.0
   kahan_c = 0.0
   do i_s = 1, n_seg
      do k = 1, nxys(i_s)
         kahan_y = slip(k, i_s)*cniu(k, i_s)-kahan_c 
         kahan_t = dsum+kahan_y
         kahan_c = (kahan_t-dsum)-kahan_y
         dsum = kahan_t
!         dsum = dsum+slip(k, i_s)*cniu(k, i_s)    ! we may need to increase numerical precision
      end do
   end do

   call define_slip_field(slip, rake)
   call tlap(rupt_time, err_time)
!
!  An iteration of the simulated annealing algorithm, over each subfault.
!
   ll = 0
!       begin to perturb       
!
   do k = 1, nnn
      nsub(k) = k
   end do

   do k = 1, nnn-1
      nran = k
      do while (nran .eq. k .or. nran .gt. nnn)
         x = ran1()
         nran = int(x*(nnn-k)+k+1)
      end do
      nbb = nsub(nran)
      nn = nsub(k)
      nsub(k) = nbb
      nsub(nran) = nn
   end do

   do k = 1, nnn
      ll = nsub(k)
      if (ll .gt. nnn) stop
      n_total = 0
      do i_ss = 1, n_seg
         n_total = nxys(i_ss)+n_total
         if (ll .le. n_total) then
            i_s = i_ss
            ll_s = ll
            exit
         end if
      end do
      do i_ss = 1, i_s-1
         ll_s = ll_s-nxys(i_ss)
      end do
      isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
      isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
      rake2 = rake(ll_s, i_s)*dpi
      a = sin(rake2)*slip(ll_s, i_s)
      b = cos(rake2)*slip(ll_s, i_s)
      a_sub = rake(ll_s, i_s)
!
!  make up unchange graph
!
      do ir = 1, n_chan
         df = df_ref/dt_channel(ir)
         omega = -twopi*df*rupt_time(ll_s, i_s)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, jfmax
!            omega = -twopi_0*df*(i-1)*rupt_time(ll_s, i_s)
!            z = cmplx(cos(omega), sin(omega))
            green_subf = (a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll))* &
            &    source(i, ir, isl, isr)*z
            forward_real(i, ir) = forward_real(i, ir)-real(green_subf)
            forward_imag(i, ir) = forward_imag(i, ir)-aimag(green_subf)
            z = z*z1    ! we may need to increase numerical precision
         end do
      end do
      dsum = dsum-slip(ll_s, i_s)*cniu(ll_s, i_s)
      nn_sub = ll
!  
      n_accept = 0
      npb = np(4*(ll-1)+1)
      if (npb .lt. 2) exit
!
!  slip extreme values
!
      npb = np(4*(ll-1)+1)
      dpb = dp(4*(ll-1)+1)
      slip_beg = beg(4*(ll-1)+1)
      slip_max = (npb-1)*dpb
      slip_end = slip_beg+slip_max
!
!  rake extreme values
!  
      npb = np(4*(ll-1)+2)
      dpb = dp(4*(ll-1)+2)
      angle_beg = beg(4*(ll-1)+2)
      angle_max = (npb-1)*dpb
      angle_end = angle_beg+angle_max
! 
!  rupture time extreme values.
!
      rupt_beg = time_min(ll_s, i_s)
      rupt_end = time_max(ll_s, i_s)
      rupt_max = rupt_end-rupt_beg
!
!  rise time parameters extreme values
!  
      rise_time_beg = ta0
      rise_time_end = ta0+(msou-1)*dta
      rise_time_max = (msou-1)*dta
      do i_move = 1, max_move
!
!       Save values before the perturbation
!
         t_save = rupt_time(ll_s, i_s)
         d_save = slip(ll_s, i_s)
         a_save = rake(ll_s, i_s)
         l_save = tl(ll_s, i_s)
         r_save = tr(ll_s, i_s)
!
!  Perturb the slip
!
         call cauchy(t, c)                           
         duse = d_save+c*slip_max
         duse = min(duse, slip_end)
         duse = max(duse, slip_beg)
!
!  Perturb the rake
!
         call cauchy(t, c)                           
         ause = a_save+c*angle_max
         ause = min(ause, angle_end)
         ause = max(ause, angle_beg)
! 
!  Perturb the rupture time.
!
         if (rupt_max .gt. min_dt) then
            call cauchy(t, c)                           
            vuse = min_dt*int((t_save+c*rupt_max)/min_dt+0.5)     
            vuse = min(vuse, rupt_end)
            vuse = max(vuse, rupt_beg)
         else
            vuse = rupt_beg
         end if             
!
!  Perturb rise time parameters
!  
         call cauchy(t, c)                           
         isl = int((l_save+c*rise_time_max-ta0)/dta+0.5)+1
         isl = min0(isl, msou)
         isl = max0(isl, 1)
         call cauchy(t, c)                           
         isr = int((r_save+c*rise_time_max-ta0)/dta+0.5)+1
         isr = min0(isr, msou)
         isr = max0(isr, 1)
         
         rake2 = ause*dpi
         a = duse*sin(rake2)
         b = duse*cos(rake2)
         ersum2 = 0.d0
!$omp parallel & 
!$omp& default(shared) &
!$omp& private(ir, df, i, omega, z, z1, green_subf, cr, cz, forward2, ex)
!$omp do schedule(guided) reduction(+:ersum2)
         do ir = 1, n_chan
            df = df_ref/dt_channel(ir)
            omega = -twopi*df*vuse
            z1 = cmplx(cos(omega), sin(omega), double) 
            z = cmplx(1.d0, 0.d0, double)
            do i = 1, jfmax
!               omega = -twopi_0*df*(i-1)*vuse
!               z = cmplx(cos(omega), sin(omega)) 
               green_subf = (a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll))* &
               &    source(i, ir, isl, isr)*z
               cr(i) = real(green_subf)+forward_real(i, ir)
               cz(i) = aimag(green_subf)+forward_imag(i, ir)
               z = z*z1    ! we may need to increase numerical precision
            end do
            call wavelet_syn(cr, cz, forward2)
            call misfit(ir, forward2, ex)     
            ersum2 = ersum2 + ex    ! we may need to increase numerical precision
         end do
!$omp end do
!$omp end parallel
         dsum = dsum+duse*cniu(ll_s, i_s)
         cmoment = dsum*area
         derr = (cmoment/cm_point)
         amp = 1.0
         call modify_slip_field(nn_sub, duse, ause)
         call lap(err_s)
         rupt_time(ll_s, i_s) = vuse
         call tlap(rupt_time, err_time)

         erm = ersum2+derr*coef_moment+amp*err_s*coef_slip
         erm = erm+coef_time*err_time+coef_st*err_static
         dsum = dsum-duse*cniu(ll_s, i_s)
         de = erm-emin
!  
!  Now, we update the kinematic model.
!  
         rand = ran1()
         aux = exp(-de/t)
         if (aux .gt. rand) then
            emin = erm
            slip(ll_s, i_s) = duse
            rake(ll_s, i_s) = ause
            rupt_time(ll_s, i_s) = vuse
            tl(ll_s, i_s) = (isl-1)*dta+ta0
            tr(ll_s, i_s) = (isr-1)*dta+ta0
            n_accept = n_accept+1
         else
            rupt_time(ll_s, i_s) = t_save
            tl(ll_s, i_s) = l_save
            tr(ll_s, i_s) = r_save
            slip(ll_s, i_s) = d_save
            rake(ll_s, i_s) = a_save
         end if
         ermin = min(ermin, erm)
         if (n_accept .gt. accept_max) exit
      end do
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  finish the perturbation of subevent (i_s, ll_s)
!
      rake2 = rake(ll_s, i_s)*dpi
      a = sin(rake2)*slip(ll_s, i_s)
      b = cos(rake2)*slip(ll_s, i_s)
      isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
      isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
      dsum = dsum+slip(ll_s, i_s)*cniu(ll_s, i_s)
      d_sub = slip(ll_s, i_s)
      a_sub = rake(ll_s, i_s)
      nn_sub = ll
      call modify_slip_field(nn_sub, d_sub, a_sub)
      call tlap(rupt_time, err_time)
      do ir = 1, n_chan
         df = df_ref/dt_channel(ir)
         omega = -twopi*df*rupt_time(ll_s, i_s)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, jfmax
!            omega = -twopi_0*df*(i-1)*rupt_time(ll_s, i_s)
!            z = cmplx(cos(omega), sin(omega))   
            green_subf = (a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll)) &
            &   *source(i, ir, isl, isr)*z
            forward_real(i, ir) = forward_real(i, ir)+real(green_subf)
            forward_imag(i, ir) = forward_imag(i, ir)+aimag(green_subf)
            z = z*z1
         end do
      end do

   end do
   write(*,*) ermin
   end subroutine annealing_iter3 

   
   subroutine annealing_iter4(slip, rake, rupt_time, tl, tr, er, t)
   !!
   !!  We compute an iteration of the annealing method, adding static GNSS data.
   !!  Parallelized version.
   !!
   use omp_lib
   implicit none
   integer isl, isr, nn_sub, nsub(nnxy_m), n_accept, &
   & nbb, i, k, npb, nn, nran, ll_s, i_s, i_ss, ir, ll, iys, &
   & ixs, i_move, n_total
   real slip(nnxy, max_seg), rake(nnxy, max_seg), rupt_time(nnxy, max_seg), &
   & tr(nnxy, max_seg), tl(nnxy, max_seg), er, t, &
   & forward_real(npth, nnsta), forward_imag(npth, nnsta), duse, ause, vuse, &
   & de, rand, c, aux, dpb, amp, derr, erm, err_static, &
   & cmoment, d_sub, a_sub, err_s, a, b, kahan_y, kahan_c, kahan_t, &
   & err_time, t_save, a_save, d_save, x, dsum, &
   & l_save, r_save, cr(inptd), cz(inptd), forward2(inptd), &
   & slip_beg, slip_max, slip_end, angle_beg, angle_end, angle_max, &
   & rupt_beg, rupt_end, rupt_max, rise_time_beg, rise_time_end, rise_time_max
   real*8 :: omega, ersum2
   real :: df, df_ref, rake2, ex
   complex :: green_subf
   complex*16 :: z, z1, forward(npth), z0
!
   z0 = cmplx(0.d0, 0.d0, double)
   erm = 0.0
!
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Here, we compute the value of the objective function, 
!  using the input kinematic model.
!
   ll = 0
   df_ref = 1.0/(2.0**lnpt)

   do ir = 1, n_chan
      df = df_ref/dt_channel(ir)
      forward(:) = z0
      ll = 0
      do i_s = 1, n_seg
         do iys = 1, nys_sub(i_s)
            do ixs = 1, nxs_sub(i_s)
               ll = ll+1
               ll_s = (iys-1)*nxs_sub(i_s)+ixs
               rake2 = rake(ll_s, i_s)*dpi
               a = sin(rake2)*slip(ll_s, i_s)
               b = cos(rake2)*slip(ll_s, i_s)
               isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
               isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
               omega = -twopi*df*rupt_time(ll_s, i_s)
               z1 = cmplx(cos(omega), sin(omega), double)
               z = cmplx(1.d0, 0.d0, double)
               do i = 1, jfmax
                  forward(i) = forward(i) &
                  & +(a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll))*source(i, ir, isl, isr)*z
                  z = z*z1    ! we may need to increase numerical precision
               end do
            end do
         end do
      end do

      do i = 1, npth
         forward_real(i, ir) = real(forward(i))
         forward_imag(i, ir) = aimag(forward(i))
      end do
   end do
   call static_synthetic(slip, rake, nxys, err_static)

   dsum = 0.0
   kahan_y = 0.0
   kahan_t = 0.0
   kahan_c = 0.0
   do i_s = 1, n_seg
      do k = 1, nxys(i_s)
         kahan_y = slip(k, i_s)*cniu(k, i_s)-kahan_c 
         kahan_t = dsum+kahan_y
         kahan_c = (kahan_t-dsum)-kahan_y
         dsum = kahan_t
!         dsum = dsum+slip(k, i_s)*cniu(k, i_s)    ! we may need to increase numerical precision
      end do
   end do

   call define_slip_field(slip, rake)
   call tlap(rupt_time, err_time)
!
!  An iteration of the simulated annealing algorithm, over each subfault.
!
   ll = 0
!       begin to perturb       
!
   do k = 1, nnn
      nsub(k) = k
   end do

   do k = 1, nnn-1
      nran = k
      do while (nran .eq. k .or. nran .gt. nnn)
         x = ran1()
         nran = int(x*(nnn-k)+k+1)
      end do
      nbb = nsub(nran)
      nn = nsub(k)
      nsub(k) = nbb
      nsub(nran) = nn
   end do

   do k = 1, nnn
      ll = nsub(k)
      if (ll .gt. nnn) stop
      n_total = 0
      do i_ss = 1, n_seg
         n_total = nxys(i_ss)+n_total
         if (ll .le. n_total) then
            i_s = i_ss
            ll_s = ll
            exit
         end if
      end do
      do i_ss = 1, i_s-1
         ll_s = ll_s-nxys(i_ss)
      end do
      isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
      isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
      rake2 = rake(ll_s, i_s)*dpi
      a = sin(rake2)*slip(ll_s, i_s)
      b = cos(rake2)*slip(ll_s, i_s)
      a_sub = rake(ll_s, i_s)
!
!  make up unchange graph
!
      do ir = 1, n_chan
         df = df_ref/dt_channel(ir)
         omega = -twopi*df*rupt_time(ll_s, i_s)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, jfmax
!            omega = -twopi_0*df*(i-1)*rupt_time(ll_s, i_s)
!            z = cmplx(cos(omega), sin(omega))
            green_subf = (a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll))* &
            &    source(i, ir, isl, isr)*z
            forward_real(i, ir) = forward_real(i, ir)-real(green_subf)
            forward_imag(i, ir) = forward_imag(i, ir)-aimag(green_subf)
            z = z*z1    ! we may need to increase numerical precision
         end do
      end do
      call static_remove_subfault(slip, rake, i_s, ll_s)
      dsum = dsum-slip(ll_s, i_s)*cniu(ll_s, i_s)
      nn_sub = ll
!  
      n_accept = 0
      npb = np(4*(ll-1)+1)
      if (npb .lt. 2) exit
!
!  slip extreme values
!
      npb = np(4*(ll-1)+1)
      dpb = dp(4*(ll-1)+1)
      slip_beg = beg(4*(ll-1)+1)
      slip_max = (npb-1)*dpb
      slip_end = slip_beg+slip_max
!
!  rake extreme values
!  
      npb = np(4*(ll-1)+2)
      dpb = dp(4*(ll-1)+2)
      angle_beg = beg(4*(ll-1)+2)
      angle_max = (npb-1)*dpb
      angle_end = angle_beg+angle_max
! 
!  rupture time extreme values.
!
      rupt_beg = time_min(ll_s, i_s)
      rupt_end = time_max(ll_s, i_s)
      rupt_max = rupt_end-rupt_beg
!
!  rise time parameters extreme values
!  
      rise_time_beg = ta0
      rise_time_end = ta0+(msou-1)*dta
      rise_time_max = (msou-1)*dta
      do i_move = 1, max_move
!
!       Save values before the perturbation
!
         t_save = rupt_time(ll_s, i_s)
         d_save = slip(ll_s, i_s)
         a_save = rake(ll_s, i_s)
         l_save = tl(ll_s, i_s)
         r_save = tr(ll_s, i_s)
!
!  Perturb the slip
!
         call cauchy(t, c)                           
         duse = d_save+c*slip_max
         duse = min(duse, slip_end)
         duse = max(duse, slip_beg)
!
!  Perturb the rake
!
         call cauchy(t, c)                           
         ause = a_save+c*angle_max
         ause = min(ause, angle_end)
         ause = max(ause, angle_beg)
! 
!  Perturb the rupture time.
!
         if (rupt_max .gt. min_dt) then
            call cauchy(t, c)                           
            vuse = min_dt*int((t_save+c*rupt_max)/min_dt+0.5)     
            vuse = min(vuse, rupt_end)
            vuse = max(vuse, rupt_beg)
         else
            vuse = rupt_beg
         end if             
!
!  Perturb rise time parameters
!  
         call cauchy(t, c)                           
         isl = int((l_save+c*rise_time_max-ta0)/dta+0.5)+1
         isl = min0(isl, msou)
         isl = max0(isl, 1)
         call cauchy(t, c)                           
         isr = int((r_save+c*rise_time_max-ta0)/dta+0.5)+1
         isr = min0(isr, msou)
         isr = max0(isr, 1)
         
         rake2 = ause*dpi
         a = duse*sin(rake2)
         b = duse*cos(rake2)
         ersum2 = 0.d0
!$omp parallel & 
!$omp& default(shared) &
!$omp& private(ir, df, i, omega, z, z1, green_subf, cr, cz, forward2, ex)
!$omp do schedule(guided) reduction(+:ersum2)
         do ir = 1, n_chan
            df = df_ref/dt_channel(ir)
            omega = -twopi*df*vuse
            z1 = cmplx(cos(omega), sin(omega), double)
            z = cmplx(1.d0, 0.d0, double)
            do i = 1, jfmax
               green_subf = (a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll))* &
               &    source(i, ir, isl, isr)*z
               cr(i) = real(green_subf)+forward_real(i, ir)
               cz(i) = aimag(green_subf)+forward_imag(i, ir)
               z = z*z1    ! we may need to increase numerical precision
            end do
            call wavelet_syn(cr, cz, forward2)
            call misfit(ir, forward2, ex)     
            ersum2 = ersum2 + ex    ! we may need to increase numerical precision
         end do
!$omp end do
!$omp end parallel
         dsum = dsum+duse*cniu(ll_s, i_s)
         cmoment = dsum*area
         derr = (cmoment/cm_point)
         amp = 1.0
         call static_modify_subfault(duse, ause, i_s, ll_s, err_static)
         call modify_slip_field(nn_sub, duse, ause)
         call lap(err_s)
         rupt_time(ll_s, i_s) = vuse
         call tlap(rupt_time, err_time)

         erm = ersum2+derr*coef_moment+amp*err_s*coef_slip
         erm = erm+coef_time*err_time+coef_st*err_static
         dsum = dsum-duse*cniu(ll_s, i_s)
         de = erm-emin
!  
!  Now, we update the kinematic model.
!  
         rand = ran1()
         aux = exp(-de/t)
         if (aux .gt. rand) then
            emin = erm
            slip(ll_s, i_s) = duse
            rake(ll_s, i_s) = ause
            rupt_time(ll_s, i_s) = vuse
            tl(ll_s, i_s) = (isl-1)*dta+ta0
            tr(ll_s, i_s) = (isr-1)*dta+ta0
            n_accept = n_accept+1
         else
            rupt_time(ll_s, i_s) = t_save
            tl(ll_s, i_s) = l_save
            tr(ll_s, i_s) = r_save
            slip(ll_s, i_s) = d_save
            rake(ll_s, i_s) = a_save
         end if
         ermin = min(ermin, erm)
         if (n_accept .gt. accept_max) exit
      end do
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  finish the perturbation of subevent (i_s, ll_s)
!
      rake2 = rake(ll_s, i_s)*dpi
      a = sin(rake2)*slip(ll_s, i_s)
      b = cos(rake2)*slip(ll_s, i_s)
      isl = int((tl(ll_s, i_s)-ta0)/dta+0.5)+1
      isr = int((tr(ll_s, i_s)-ta0)/dta+0.5)+1
      dsum = dsum+slip(ll_s, i_s)*cniu(ll_s, i_s)
      d_sub = slip(ll_s, i_s)
      a_sub = rake(ll_s, i_s)
      nn_sub = ll
      call modify_slip_field(nn_sub, d_sub, a_sub)
      call tlap(rupt_time, err_time)
      do ir = 1, n_chan
         df = df_ref/dt_channel(ir)
         omega = -twopi*df*rupt_time(ll_s, i_s)
         z1 = cmplx(cos(omega), sin(omega), double)   
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, jfmax
!            omega = -twopi_0*df*(i-1)*rupt_time(ll_s, i_s)
!            z = cmplx(cos(omega), sin(omega))   
            green_subf = (a*green_dip(i, ir, ll)+b*green_stk(i, ir, ll)) &
            &   *source(i, ir, isl, isr)*z
            forward_real(i, ir) = forward_real(i, ir)+real(green_subf)
            forward_imag(i, ir) = forward_imag(i, ir)+aimag(green_subf)
            z = z*z1
         end do
      end do
      call static_add_subfault(slip, rake, i_s, ll_s, err_static)

   end do
   write(*,*) ermin
   end subroutine annealing_iter4


end module annealing 
