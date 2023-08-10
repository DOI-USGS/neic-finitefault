module annealing


   use constants, only : pi, max_stations, dpi, twopi, max_subf, max_seg, wave_pts2, &
            &   max_subfaults2, wave_pts, max_subfaults
   use retrieve_gf, only : green_dip, green_stk
   use wavelets, only : wavelet_syn
   use wavelet_param, only : get_data_param 
   use rise_time, only : source
   use get_stations_data, only : get_properties, count_wavelets
   use random_gen, only : ran1, cauchy
   use misfit_eval, only : misfit_channel
   use modelling_inputs, only : smooth_moment, smooth_slip, smooth_time, io_re, moment_input, emin0, &
            &   t_latest
   use regularization, only : slip_laplace, time_laplace, define_slip_field, modify_slip_field
   use static_data, only : static_synthetic, static_remove_subfault, &
                       &   static_modify_subfault, static_add_subfault
   use insar_data, only : insar_synthetic, insar_remove_subfault, &
                       &   insar_modify_subfault, insar_add_subfault, &
                       &   insar_remove_ramp, insar_modify_ramp, &
                       &   insar_add_ramp, ramp_length
   use omp_lib
   implicit none
   real :: coef_moment, coef_slip, coef_gps, coef_insar, coef_time
   real :: current_value, min_value, min_dt, area
   real :: insar_misfit0
   integer :: subfaults_segment(max_seg)
   integer, parameter :: double = kind(1.d0)
   integer, private :: threads
   integer, parameter, private :: max_move=50, accept_max=5
   integer :: segments, msou, nxs_sub(max_seg), nys_sub(max_seg), subfaults
   real :: beg(max_subfaults2), dp(max_subfaults2)
   integer :: np(max_subfaults2)
   real :: time_min(max_subfaults), time_max(max_subfaults)
   real :: shear(max_subfaults), ta0, dta, dxs, dys
   integer :: max_freq, lnpt, channels
   real :: dt_channel(max_stations)
   logical :: segment_in_event(max_seg, 10), subfault_in_event(max_subfaults, 10)
   real*8, allocatable :: forward_real2(:, :), forward_imag2(:, :)
   real, allocatable :: forward_real(:, :), forward_imag(:, :)


contains


   subroutine n_threads(auto)
!
!  Args:
!  auto: True if automatic inversion, False otherwise
!
   implicit none
   logical :: auto
   threads = 4
   if (auto) threads = 3
   call omp_set_num_threads(threads)
   end subroutine n_threads


   subroutine annealing_set_fault_parameters()
   use model_parameters, only : get_rise_time, get_shear, get_segments, &
   &  get_subfaults, get_space, get_events_segments
   implicit none
   real :: dip(max_seg), strike(max_seg), delay_seg(max_seg)
   integer :: cum_subfaults(max_seg), nx_p, ny_p
   real :: v_min, v_max, v_ref, time_ref(max_subfaults2)
   call get_rise_time(ta0, dta, msou)
   call get_shear(shear)
   call get_segments(nxs_sub, nys_sub, dip, strike, delay_seg, segments, subfaults, cum_subfaults)
   call get_subfaults(dxs, dys, nx_p, ny_p, v_min, v_max, v_ref)
   call get_space(time_min, time_max, time_ref, beg, dp, np)
   call get_events_segments(segment_in_event, subfault_in_event)
   end subroutine annealing_set_fault_parameters


   subroutine annealing_set_data_properties()
   implicit none
   integer :: jmin, jmax, nlen
   character(len=15) :: sta_name(max_stations)
   character(len=3) :: component(max_stations)

   call get_properties(sta_name, component, dt_channel, channels)
   call get_data_param(lnpt, jmin, jmax, nlen, max_freq)
   end subroutine annealing_set_data_properties


   subroutine allocate_forward()
   implicit none
   allocate(forward_real2(wave_pts, max_stations))
   allocate(forward_imag2(wave_pts, max_stations))
   allocate(forward_real(wave_pts, max_stations))
   allocate(forward_imag(wave_pts, max_stations))
   end subroutine allocate_forward


   subroutine deallocate_forward()
   implicit none
   deallocate(forward_real2)
   deallocate(forward_imag2)
   deallocate(forward_real)
   deallocate(forward_imag)
   end subroutine deallocate_forward 
   
   
   subroutine initial_model(slip, rake, rupt_time, t_rise, t_fall)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  t_rise: array with model risetime values for all subfaults
!  t_fall: array with model falltime values for all subfaults
!
   use model_parameters, only : slip0, rake0, rupt_time0, t_rise0, t_fall0
   implicit none
   real :: slip(:), rake(:), rupt_time(:), t_rise(:), t_fall(:)
   real :: x(max_subfaults2)
   integer :: nxy, i, segment, k, npa, subfault, etc
   
   npa = 0
   do segment = 1, segments
      npa = npa + 4*nxs_sub(segment)*nys_sub(segment)
   end do
   do k = 1, npa
      x(k) = beg(k)+int(ran1()*np(k)-1.0)*dp(k)
   end do
   if (io_re .eq. 0) then
      k = 0
      subfault = 0
      do segment = 1, segments
         nxy = nxs_sub(segment)*nys_sub(segment)
         subfaults_segment(segment) = nxy
         do i = 1, nxy
            subfault = subfault + 1
            k = k + 1
            slip(subfault) = x(k)
            k = k + 1
            rake(subfault) = x(k)
            k = k + 2
            t_rise(subfault) = (ta0+int(ran1()*msou)*dta)
            t_fall(subfault) = (ta0+int(ran1()*msou)*dta)
            rupt_time(subfault) = ran1()* &
     &  (time_max(subfault)-time_min(subfault))+time_min(subfault)
         end do
      end do 
   else
      slip(:) = slip0(:)
      rake(:) = rake0(:)
      rupt_time(:) = rupt_time0(:)
      t_rise(:) = t_rise0(:)
      t_fall(:) = t_fall0(:)
   end if
   end subroutine initial_model


   subroutine print_summary(slip, rake, rupt_time, t_rise, t_fall, static, insar, &
      &   get_coeff, ramp)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  t_rise: array with model risetime values for all subfaults
!  t_fall: array with model falltime values for all subfaults
!  static: True if static GPS data used in modelling, False otherwise
!  insar: True if insar data used in modelling, False otherwise
!  get_coeff: get regularization coefficients if and only if this is True
!  ramp: Value of insar ramp, optional
!
   implicit none
   real*8, optional :: ramp(:)
   real :: slip(:), rake(:), rupt_time(:)
   real :: t_rise(:), t_fall(:)
   real amp, moment, moment_reg, dt, value1, slip_reg, gps_misfit, insar_misfit, &
      & time_reg, forward_real3(wave_pts2, max_stations), a, b, &
      & forward_imag3(wave_pts2, max_stations), real1(wave_pts2), imag1(wave_pts2), coeffs_syn(wave_pts2)
   real :: rake2, delta_freq, delta_freq0, moment0, kahan_y, kahan_t, kahan_c
   real*8 :: omega, misfit2, ex
   integer :: i, segment, channel, isl, isr, jf, k, subfault, used_data
   complex*16 :: z0, forward(wave_pts), z, z1
   logical :: static, get_coeff, insar
   character(len=15) :: sta_name(max_stations)
   character(len=3) :: component(max_stations)

   used_data = 0
   call count_wavelets(used_data)
   z0 = cmplx(0.d0, 0.d0, double)
   min_dt = 10
   gps_misfit = 0.0
   insar_misfit = 0.0
   do i = 1, max_stations
      if (dt_channel(i) .gt. 1.e-4) min_dt = min(min_dt, dt_channel(i))
   end do
   jf = 2**(lnpt-1)+1
!
! Compute synthetics given current fault model
!
   area = dxs*dys*(1.e+10)
   delta_freq0 = 1.0/(2.0**lnpt)
   misfit2 = 0.d0
   do channel = 1, channels
      delta_freq = delta_freq0/dt_channel(channel)
      dt = dt_channel(channel)
      do i = 1, wave_pts
         forward(i) = z0
      end do
      do subfault = 1, subfaults
         rake2 = rake(subfault)*dpi
         a = sin(rake2)*slip(subfault)
         b = cos(rake2)*slip(subfault)
         isl = int((t_rise(subfault)-ta0)/dta+0.5)+1
         isr = int((t_fall(subfault)-ta0)/dta+0.5)+1
         omega = -twopi*delta_freq*rupt_time(subfault)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, max_freq
            forward(i) = forward(i) &
         &  +(a*green_dip(i, channel, subfault)+b*green_stk(i, channel, subfault)) &
         &  *source(i, channel, isl, isr)*z
            z = z*z1
         end do
      end do

      do i = 1, wave_pts
         forward_real3(i, channel) = real(forward(i))
         forward_imag3(i, channel) = aimag(forward(i))
         real1(i) = forward_real3(i, channel)
         imag1(i) = forward_imag3(i, channel)
      end do
      call wavelet_syn(real1, imag1, coeffs_syn)
      call misfit_channel(channel, coeffs_syn, ex)
      misfit2 = misfit2 + ex
   end do

   amp = 1.0
   moment0 = 0.0
   kahan_y = 0.0
   kahan_t = 0.0
   kahan_c = 0.0
   do subfault = 1, subfaults
      kahan_y = slip(subfault)*shear(subfault)-kahan_c 
      kahan_t = moment0+kahan_y
      kahan_c = (kahan_t-moment0)-kahan_y
      moment0 = kahan_t
!         moment0 = moment0+slip(k, segment)*shear(k, segment)    ! we may need to increase numerical precision
   end do

   current_value = emin0
   moment = moment0*area
   moment_reg = (moment/moment_input)-1
   if(moment_reg.ge. 0.10)then
      moment_reg = sqrt(5*moment_reg+0.5)
   elseif((-0.1 .lt. moment_reg) .and. (moment_reg .lt. 0.1))then
      moment_reg = (10*moment_reg)**4
   else
      moment_reg = sqrt(-5*moment_reg+0.5)
   endif
!   moment_reg = (moment/moment_input)
   call define_slip_field(slip, rake)
   call slip_laplace(slip_reg)
   call time_laplace(rupt_time, time_reg)
   if (time_reg .lt. 1.0e-9) time_reg = 1.0
 
   coef_gps = 0.0
   coef_insar = 0.0
   insar_misfit0 = 0.0
   if (static) call static_synthetic(slip, rake, gps_misfit)
   if (insar) then
      if (present(ramp)) then
         call insar_synthetic(slip, rake, insar_misfit, ramp)
      else
         call insar_synthetic(slip, rake, insar_misfit)
      endif
      insar_misfit0 = insar_misfit
   endif

   if (get_coeff) then
      coef_moment=smooth_moment*(misfit2 - current_value)
      coef_moment = min(coef_moment, 1.0)
      coef_slip = smooth_slip*misfit2/(slip_reg*amp)
      coef_slip = min(0.003, coef_slip)
      coef_time = smooth_time*misfit2/(time_reg*amp)
      if (static) then
         coef_gps = misfit2/(gps_misfit*amp)
      endif
      if (insar) then
         coef_insar = misfit2/(insar_misfit*amp)
      endif
   endif

   value1 = misfit2+coef_moment*moment_reg+coef_slip*slip_reg*amp
   value1 = value1+coef_time*time_reg
   value1 = value1+coef_gps*gps_misfit+coef_insar*insar_misfit
   do segment = 1, segments
      subfaults_segment(segment) = nys_sub(segment)*nxs_sub(segment)
   end do
   write(*,'()')
   write(*,*)'averaged misfit error', misfit2
   write(*,*)'moment error', moment_reg
   write(*,*)'slip smoothness penalization', slip_reg
   write(*,*)'time smoothness penalization', time_reg
   if (static) write(*,*)'static data penalization', gps_misfit
   if (insar) write(*,*)'insar data penalization', insar_misfit
   write(*,*)'total moment of the inversion', moment
   write(*,'(/A, F10.7)')'moment error coefficient', coef_moment
   write(*,*)'slip smoothness penalization coefficient', coef_slip
   write(*,*)'time smoothness penalization coefficient', coef_time
   if (static) write(*,*)'static data penalization coefficient', coef_gps
   if (insar) write(*,*)'insar data penalization coefficient', coef_insar
   write(*,'(/A, I4)')'Amount of variables: ', 5 * subfaults
   write(*,*)'Amount of data values: ', used_data
   current_value = value1
   min_value = value1
   open(12,file='modelling_summary.txt')
   write(12,'(/A/)')'Modelling Report'
   write(12,*)'averaged misfit error', misfit2
   write(12,*)'moment error', moment_reg
   write(12,*)'slip smoothness penalization', slip_reg
   write(12,*)'time smoothness penalization', time_reg
   if (static) write(12,*)'static data penalization', gps_misfit
   if (insar) write(12,*)'insar data penalization', insar_misfit
   write(12,*)'objective function value', value1
   write(12,*)'total moment of the inversion', moment
   write(12,'(/A, F10.7)')'moment error coefficient', coef_moment
   write(12,*)'slip smoothness penalization coefficient', coef_slip
   write(12,*)'time smoothness penalization coefficient', coef_time
   if (static) write(12,*)'static data penalization coefficient', coef_gps
   if (insar) write(12,*)'insar data penalization coefficient', coef_insar
   write(12,'(/A, I4)')'Amount of variables: ', 5 * subfaults
   write(12,*)'Amount of data values: ', used_data
   close(12)
   end subroutine print_summary
   
   
   subroutine print_summary2(slip, rake, rupt_time, t_rise, t_fall, static, insar, &
      &   get_coeff, ramp)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  t_rise: array with model risetime values for all subfaults
!  t_fall: array with model falltime values for all subfaults
!  static: True if static GPS data used in modelling, False otherwise
!  insar: True if insar data used in modelling, False otherwise
!  get_coeff: get regularization coefficients if and only if this is True
!  ramp: Value of insar ramp, optional
!
   use modelling_inputs, only : events, moment_event
   implicit none
   real*8, optional :: ramp(:)
   real :: slip(:), rake(:), rupt_time(:)
   real :: t_rise(:), t_fall(:)
   real amp, moment, moment_reg(10), moment_reg2, dt, value1, slip_reg, gps_misfit, insar_misfit, &
      & time_reg, forward_real3(wave_pts2, max_stations), a, b, &
      & forward_imag3(wave_pts2, max_stations), real1(wave_pts2), imag1(wave_pts2), coeffs_syn(wave_pts2)
   real :: rake2, delta_freq, delta_freq0, moment0(10), kahan_y, kahan_t, kahan_c
   real*8 :: omega, misfit2, ex
   integer :: i, segment, channel, isl, isr, jf, k, subfault, used_data
   complex*16 :: z0, forward(wave_pts), z, z1
   logical :: static, get_coeff, insar
   character(len=15) :: sta_name(max_stations)
   character(len=3) :: component(max_stations)

   used_data = 0
   call count_wavelets(used_data)
   z0 = cmplx(0.d0, 0.d0, double)
   min_dt = 10
   gps_misfit = 0.0
   insar_misfit = 0.0
   do i = 1, max_stations
      if (dt_channel(i) .gt. 1.e-4) min_dt = min(min_dt, dt_channel(i))
   end do
   jf = 2**(lnpt-1)+1
!
! Compute synthetics given current fault model
!
   area = dxs*dys*(1.e+10)
   delta_freq0 = 1.0/(2.0**lnpt)
   misfit2 = 0.d0
   do channel = 1, channels
      delta_freq = delta_freq0/dt_channel(channel)
      dt = dt_channel(channel)
      do i = 1, wave_pts
         forward(i) = z0
      end do
      do subfault = 1, subfaults
         rake2 = rake(subfault)*dpi
         a = sin(rake2)*slip(subfault)
         b = cos(rake2)*slip(subfault)
         isl = int((t_rise(subfault)-ta0)/dta+0.5)+1
         isr = int((t_fall(subfault)-ta0)/dta+0.5)+1
         omega = -twopi*delta_freq*rupt_time(subfault)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, max_freq
            forward(i) = forward(i) &
         &  +(a*green_dip(i, channel, subfault)+b*green_stk(i, channel, subfault)) &
         &  *source(i, channel, isl, isr)*z
            z = z*z1
         end do
      end do

      do i = 1, wave_pts
         forward_real3(i, channel) = real(forward(i))
         forward_imag3(i, channel) = aimag(forward(i))
         real1(i) = forward_real3(i, channel)
         imag1(i) = forward_imag3(i, channel)
      end do
      call wavelet_syn(real1, imag1, coeffs_syn)
      call misfit_channel(channel, coeffs_syn, ex)
      misfit2 = misfit2 + ex
   end do

   amp = 1.0
   do k=1, events
      moment0(k) = 0.0
      kahan_y = 0.0
      kahan_t = 0.0
      kahan_c = 0.0
      do subfault = 1, subfaults
         if (subfault_in_event(subfault, k)) then
            kahan_y = slip(subfault)*shear(subfault)-kahan_c 
            kahan_t = moment0(k)+kahan_y
            kahan_c = (kahan_t-moment0(k))-kahan_y
            moment0(k) = kahan_t
!         moment0 = moment0+slip(k, segment)*shear(k, segment)    ! we may need to increase numerical precision
         endif
      enddo
      moment = moment0(k)*area
      moment_reg(k) = (moment/moment_event(k))-1
      if(abs(moment_reg(k)) .ge. 0.10)then
         moment_reg(k) = sqrt(5*abs(moment_reg(k))+0.5)
      else
         moment_reg(k) = (10*abs(moment_reg(k)))**4
      endif
   end do
   moment_reg2 = 0.0
   moment = 0.0
   do k=1, events
      moment_reg2 = moment_reg2 + moment_reg(k)
      moment = moment + moment0(k)*area
   enddo

   current_value = emin0
!   moment_reg = (moment/moment_input)
   call define_slip_field(slip, rake)
   call slip_laplace(slip_reg)
   call time_laplace(rupt_time, time_reg)
   if (time_reg .lt. 1.0e-9) time_reg = 1.0
 
   coef_gps = 0.0
   coef_insar = 0.0
   insar_misfit0 = 0.0
   if (static) call static_synthetic(slip, rake, gps_misfit)
   if (insar) then
      if (present(ramp)) then
         call insar_synthetic(slip, rake, insar_misfit, ramp)
      else
         call insar_synthetic(slip, rake, insar_misfit)
      endif
      insar_misfit0 = insar_misfit
   endif

   if (get_coeff) then
      coef_moment=smooth_moment*(misfit2 - current_value)
      coef_moment = min(coef_moment, 1.0)
      coef_slip = smooth_slip*misfit2/(slip_reg*amp)
      coef_slip = min(0.003, coef_slip)
      coef_time = smooth_time*misfit2/(time_reg*amp)
      if (static) then
         coef_gps = misfit2/(gps_misfit*amp)
      endif
      if (insar) then
         coef_insar = misfit2/(insar_misfit*amp)
      endif
   endif

   value1 = misfit2+coef_moment*moment_reg2+coef_slip*slip_reg*amp
   value1 = value1+coef_time*time_reg
   value1 = value1+coef_gps*gps_misfit+coef_insar*insar_misfit
   do segment = 1, segments
      subfaults_segment(segment) = nys_sub(segment)*nxs_sub(segment)
   end do
   write(*,'()')
   write(*,*)'averaged misfit error', misfit2
   do k=1, events
      write(*,*)'moment error', moment_reg(k)
   enddo
   write(*,*)'slip smoothness penalization', slip_reg
   write(*,*)'time smoothness penalization', time_reg
   if (static) write(*,*)'static data penalization', gps_misfit
   if (insar) write(*,*)'insar data penalization', insar_misfit
   write(*,*)'total moment of the inversion', moment
   write(*,'(/A, F10.7)')'moment error coefficient', coef_moment
   write(*,*)'slip smoothness penalization coefficient', coef_slip
   write(*,*)'time smoothness penalization coefficient', coef_time
   if (static) write(*,*)'static data penalization coefficient', coef_gps
   if (insar) write(*,*)'insar data penalization coefficient', coef_insar
   write(*,'(/A, I4)')'Amount of variables: ', 5 * subfaults
   write(*,*)'Amount of data values: ', used_data
   current_value = value1
   min_value = value1
   open(12,file='modelling_summary.txt')
   write(12,'(/A/)')'Modelling Report'
   write(12,*)'averaged misfit error', misfit2
   do k=1, events
      write(12,*)'moment error', moment_reg(k)
   enddo
   write(12,*)'slip smoothness penalization', slip_reg
   write(12,*)'time smoothness penalization', time_reg
   if (static) write(12,*)'static data penalization', gps_misfit
   if (insar) write(12,*)'insar data penalization', insar_misfit
   write(12,*)'objective function value', value1
   write(12,*)'total moment of the inversion', moment
   write(12,'(/A, F10.7)')'moment error coefficient', coef_moment
   write(12,*)'slip smoothness penalization coefficient', coef_slip
   write(12,*)'time smoothness penalization coefficient', coef_time
   if (static) write(12,*)'static data penalization coefficient', coef_gps
   if (insar) write(12,*)'insar data penalization coefficient', coef_insar
   write(12,'(/A, I4)')'Amount of variables: ', 5 * subfaults
   write(12,*)'Amount of data values: ', used_data
   close(12)
   end subroutine print_summary2


   subroutine annealing_iter3(slip, rake, rupt_time, t_rise, t_fall, t)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  t_rise: array with model risetime values for all subfaults
!  t_fall: array with model falltime values for all subfaults
!  t: Current temperature of the annealing method
!
   implicit none
   real, intent(inout) :: slip(:), rake(:), rupt_time(:), t_fall(:), t_rise(:)
   real, intent(in) :: t
   integer isl, isr, n_subfault(max_subfaults), n_accept, &
   & nbb, i, k, npb, nn, nran, subfault_seg, segment, channel, subfault, &
   & n_total, j
!   real, allocatable :: forward_real(:, :), forward_imag(:, :)
   real :: duse, ause, vuse, &
   & diff, rand, c, aux, dpb, amp, moment_reg, value1, gps_misfit, &
   & moment, d_sub, a_sub, slip_reg, a, b, moment0, &
   & time_reg, t_save, a_save, d_save, x, kahan_y, kahan_t, kahan_c, &
   & l_save, r_save, real1(wave_pts2), imag1(wave_pts2), coeffs_syn(wave_pts2), &
   & slip_beg, slip_max, slip_end, angle_beg, angle_end, angle_max, &
   & rupt_beg, rupt_end, rupt_max, rise_time_beg, rise_time_end, rise_time_max
!   real*8, allocatable :: forward_real2(:, :), forward_imag2(:, :)
   real :: delta_freq0, delta_freq, rake2!, ex!, misfit2
   real*8 :: omega, misfit2, ex
   complex :: green_subf
   complex*16 :: z, z1, forward(wave_pts), z0

   z0 = cmplx(0.d0, 0.d0, double)
   value1 = 0.0
   gps_misfit = 0.0
!
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Here, we compute the value of the objective function, 
!  using the input kinematic model.
!
   subfault = 0
   delta_freq0 = 1.0/(2.0**lnpt)

   do channel = 1, channels
      delta_freq = delta_freq0/dt_channel(channel)
      forward(:) = z0
      do subfault = 1, subfaults
         rake2 = rake(subfault)*dpi
         a = sin(rake2)*slip(subfault)
         b = cos(rake2)*slip(subfault)
         isl = int((t_rise(subfault)-ta0)/dta+0.5)+1
         isr = int((t_fall(subfault)-ta0)/dta+0.5)+1
         omega = -twopi*delta_freq*rupt_time(subfault)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, max_freq
!            omega = -twopi*delta_freq*rupt_time(subfault)*(i-1)
!            z = cmplx(cos(omega), sin(omega))
            forward(i) = forward(i) &
            & +(a*green_dip(i, channel, subfault)+b*green_stk(i, channel, subfault)) &
            & *source(i, channel, isl, isr)*z
            z = z*z1    ! we may need to increase numerical precision
         end do
      end do

      do i = 1, wave_pts
         forward_real2(i, channel) = dble(forward(i))
         forward_imag2(i, channel) = aimag(forward(i))
      end do
   end do

   moment0 = 0.0
   kahan_y = 0.0
   kahan_t = 0.0
   kahan_c = 0.0
   do subfault = 1, subfaults
      kahan_y = slip(subfault)*shear(subfault)-kahan_c 
      kahan_t = moment0+kahan_y
      kahan_c = (kahan_t-moment0)-kahan_y
      moment0 = kahan_t
!         moment0 = moment0+slip(k, segment)*shear(k, segment)    ! we may need to increase numerical precision
   end do

   call define_slip_field(slip, rake)
   call time_laplace(rupt_time, time_reg)
!
!  An iteration of the simulated annealing algorithm, over each subfault.
!
!       begin to perturb       
!
   do k = 1, subfaults
      n_subfault(k) = k
   end do

   do k = 1, subfaults-1
      nran = k
      do while (nran .eq. k .or. nran .gt. subfaults)
         x = ran1()
         nran = int(x*(subfaults-k)+k+1)
      end do
      nbb = n_subfault(nran)
      nn = n_subfault(k)
      n_subfault(k) = nbb
      n_subfault(nran) = nn
   end do

   do k = 1, subfaults
      subfault = n_subfault(k)
      if (subfault .gt. subfaults) stop
      n_total = 0
      do i = 1, segments
         n_total = subfaults_segment(i)+n_total
         if (subfault .le. n_total) then
            segment = i
            subfault_seg = subfault
            exit
         end if
      end do
      do i = 1, segment-1
         subfault_seg = subfault_seg-subfaults_segment(i)
      end do
      isl = int((t_rise(subfault)-ta0)/dta+0.5)+1
      isr = int((t_fall(subfault)-ta0)/dta+0.5)+1
      rake2 = rake(subfault)*dpi
      a = sin(rake2)*slip(subfault)
      b = cos(rake2)*slip(subfault)
!
!  make up unchange graph
!
      do channel = 1, channels
         delta_freq = delta_freq0/dt_channel(channel)
         omega = -twopi*delta_freq*rupt_time(subfault)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, max_freq
!            omega = -twopi*delta_freq*(i-1)*rupt_time(subfault)
!            z = cmplx(cos(omega), sin(omega))
            green_subf = &
            & (a*green_dip(i, channel, subfault)+b*green_stk(i, channel, subfault))* &
            & source(i, channel, isl, isr)*z
!            forward_real(i, channel) = forward_real(i, channel)-real(green_subf)
!            forward_imag(i, channel) = forward_imag(i, channel)-aimag(green_subf)
            forward_real2(i, channel) = forward_real2(i, channel)-real(green_subf)
            forward_imag2(i, channel) = forward_imag2(i, channel)-aimag(green_subf)
            forward_real(i, channel) = real(forward_real2(i, channel))
            forward_imag(i, channel) = real(forward_imag2(i, channel))
            z = z*z1    ! we may need to increase numerical precision
         end do
      end do
!      moment0 = moment0-slip(subfault_seg, segment)*shear(subfault_seg, segment)
      kahan_y = -slip(subfault)*shear(subfault)-kahan_c 
      kahan_t = moment0+kahan_y
      kahan_c = (kahan_t-moment0)-kahan_y
      moment0 = kahan_t
!  
      n_accept = 0
      npb = np(4*(subfault-1)+1)
      if (npb .lt. 2) exit
!
!  slip extreme values
!
      npb = np(4*(subfault-1)+1)
      dpb = dp(4*(subfault-1)+1)
      slip_beg = beg(4*(subfault-1)+1)
      slip_max = (npb-1)*dpb
      slip_end = slip_beg+slip_max
!
!  rake extreme values
!  
      npb = np(4*(subfault-1)+2)
      dpb = dp(4*(subfault-1)+2)
      angle_beg = beg(4*(subfault-1)+2)
      angle_max = (npb-1)*dpb
      angle_end = angle_beg+angle_max
! 
!  rupture time extreme values.
!
      rupt_beg = time_min(subfault)
      rupt_end = time_max(subfault)
      rupt_max = rupt_end-rupt_beg
!
!  rise time parameters extreme values
!  
      rise_time_beg = ta0
      rise_time_end = ta0+(msou-1)*dta
      rise_time_max = (msou-1)*dta
      do i = 1, max_move
!
!       Save values before the perturbation
!
         t_save = rupt_time(subfault)
         d_save = slip(subfault)
         a_save = rake(subfault)
         l_save = t_rise(subfault)
         r_save = t_fall(subfault)
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
         vuse = rupt_beg - 0.1
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
         misfit2 = 0.d0
!$omp parallel & 
!$omp& default(shared) &
!$omp& private(channel, delta_freq, j, omega, z, z1, green_subf, real1, imag1, coeffs_syn, ex)
!$omp do schedule(static) reduction(+:misfit2)
         do channel = 1, channels
            delta_freq = delta_freq0/dt_channel(channel)
            omega = -twopi*delta_freq*vuse
            z1 = cmplx(cos(omega), sin(omega), double) 
            z = cmplx(1.d0, 0.d0, double)
            do j = 1, max_freq
!               omega = -twopi_0*delta_freq*(i-1)*vuse
!               z = cmplx(cos(omega), sin(omega)) 
               green_subf = &
               & (a*green_dip(j, channel, subfault)+b*green_stk(j, channel, subfault))* &
               & source(j, channel, isl, isr)*z
               real1(j) = real(green_subf)+forward_real(j, channel)
               imag1(j) = aimag(green_subf)+forward_imag(j, channel)
               z = z*z1    ! we may need to increase numerical precision
            end do
            call wavelet_syn(real1, imag1, coeffs_syn)
            call misfit_channel(channel, coeffs_syn, ex)     
            misfit2 = misfit2 + ex    ! we may need to increase numerical precision
         end do
!$omp end do
!$omp end parallel
         moment0 = moment0+duse*shear(subfault)
         moment = moment0*area
         moment_reg = (moment/moment_input)-1
         if(abs(moment_reg) .ge. 0.10)then
            moment_reg = sqrt(5*abs(moment_reg)+0.5)
         else
            moment_reg = (10*abs(moment_reg))**4
         endif
!         moment_reg = (moment/moment_input)
         amp = 1.0
         call modify_slip_field(subfault, duse, ause)
         call slip_laplace(slip_reg)
         rupt_time(subfault) = vuse
         call time_laplace(rupt_time, time_reg)

         value1 = misfit2+moment_reg*coef_moment+amp*slip_reg*coef_slip
         value1 = value1+coef_time*time_reg+coef_gps*gps_misfit
         moment0 = moment0-duse*shear(subfault)
         diff = value1-current_value
!  
!  Now, we update the kinematic model.
!  
         rand = ran1()
         aux = exp(-diff/t)
         if (aux .gt. rand) then
            current_value = value1
            slip(subfault) = duse
            rake(subfault) = ause
            rupt_time(subfault) = vuse
            t_rise(subfault) = (isl-1)*dta+ta0
            t_fall(subfault) = (isr-1)*dta+ta0
            n_accept = n_accept+1
         else
            rupt_time(subfault) = t_save
            t_rise(subfault) = l_save
            t_fall(subfault) = r_save
            slip(subfault) = d_save
            rake(subfault) = a_save
         end if
         min_value = min(min_value, value1)
         if (n_accept .gt. accept_max) exit
      end do
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  finish the perturbation of subevent (segment, subfault_seg)
!
      rake2 = rake(subfault)*dpi
      a = sin(rake2)*slip(subfault)
      b = cos(rake2)*slip(subfault)
      isl = int((t_rise(subfault)-ta0)/dta+0.5)+1
      isr = int((t_fall(subfault)-ta0)/dta+0.5)+1
      kahan_y = slip(subfault)*shear(subfault)-kahan_c 
      kahan_t = moment0+kahan_y
      kahan_c = (kahan_t-moment0)-kahan_y
      moment0 = kahan_t
      d_sub = slip(subfault)
      a_sub = rake(subfault)
      call modify_slip_field(subfault, d_sub, a_sub)
      do channel = 1, channels
         delta_freq = delta_freq0/dt_channel(channel)
         omega = -twopi*delta_freq*rupt_time(subfault)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, max_freq
            green_subf = (a*green_dip(i, channel, subfault)+b*green_stk(i, channel, subfault)) &
            &   *source(i, channel, isl, isr)*z
            forward_real2(i, channel) = forward_real2(i, channel)+real(green_subf)
            forward_imag2(i, channel) = forward_imag2(i, channel)+aimag(green_subf)
            forward_real(i, channel) = real(forward_real2(i, channel))
            forward_imag(i, channel) = real(forward_imag2(i, channel))
            z = z*z1
         end do
      end do

   end do
   write(*,*) min_value
   end subroutine annealing_iter3 

   
   subroutine annealing_iter4(slip, rake, rupt_time, t_rise, &
   & t_fall, t, static, insar, ramp)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  t_rise: array with model risetime values for all subfaults
!  t_fall: array with model falltime values for all subfaults
!  t: Current temperature of the annealing method
!  static: True if static GPS data used in modelling, False otherwise
!  insar: True if insar data used in modelling, False otherwise
!  ramp: Value of insar ramp, optional
!
   implicit none
   real, intent(inout) :: slip(:), rake(:), rupt_time(:), t_fall(:), t_rise(:)
   real, intent(in) :: t
   logical, intent(in) :: static, insar
   real*8, optional :: ramp(:)
   integer isl, isr, n_subfault(max_subfaults), n_accept, &
   & nbb, i, k, npb, nn, nran, subfault_seg, segment, channel, subfault, &
   & n_total, j
!   real, allocatable :: forward_real(:, :), forward_imag(:, :)
   real :: duse, ause, vuse, &
   & diff, rand, c, aux, dpb, amp, moment_reg, value1, gps_misfit, insar_misfit, &
   & moment, d_sub, a_sub, slip_reg, a, b, kahan_y, kahan_c, kahan_t, &
   & time_reg, t_save, a_save, d_save, x, moment0, &
   & l_save, r_save, real1(wave_pts2), imag1(wave_pts2), coeffs_syn(wave_pts2), &
   & slip_beg, slip_max, slip_end, angle_beg, angle_end, angle_max, &
   & rupt_beg, rupt_end, rupt_max, rise_time_beg, rise_time_end, rise_time_max
   real ramp_beg, ramp_end, ramp_max, ramp_use
   real*8 :: ramp0(36), ramp1(36)
!   real*8, allocatable :: forward_real2(:, :), forward_imag2(:, :)
   real*8 :: omega, misfit2, ex
   real :: delta_freq, delta_freq0, rake2!, ex
   complex :: green_subf
   complex*16 :: z, z1, forward(wave_pts), z0
!
   z0 = cmplx(0.d0, 0.d0, double)
   value1 = 0.0
!
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Here, we compute the value of the objective function, 
!  using the input kinematic model.
!
   subfault = 0
   delta_freq0 = 1.0/(2.0**lnpt)

   do channel = 1, channels
      delta_freq = delta_freq0/dt_channel(channel)
      forward(:) = z0
      do subfault = 1, subfaults
         rake2 = rake(subfault)*dpi
         a = sin(rake2)*slip(subfault)
         b = cos(rake2)*slip(subfault)
         isl = int((t_rise(subfault)-ta0)/dta+0.5)+1
         isr = int((t_fall(subfault)-ta0)/dta+0.5)+1
         omega = -twopi*delta_freq*rupt_time(subfault)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, max_freq
            forward(i) = forward(i) &
            & +(a*green_dip(i, channel, subfault)+b*green_stk(i, channel, subfault)) &
            & *source(i, channel, isl, isr)*z
            z = z*z1    ! we may need to increase numerical precision
         end do
      end do

      do i = 1, wave_pts
         forward_real2(i, channel) = dble(forward(i))
         forward_imag2(i, channel) = aimag(forward(i))
      end do
   end do
   if (static) call static_synthetic(slip, rake, gps_misfit)
   if (insar) then
      if (present(ramp)) then
         call insar_synthetic(slip, rake, insar_misfit, ramp)
      else
         call insar_synthetic(slip, rake, insar_misfit)
      endif
   endif
   moment0 = 0.0
   kahan_y = 0.0
   kahan_t = 0.0
   kahan_c = 0.0
   do subfault = 1, subfaults
      kahan_y = slip(subfault)*shear(subfault)-kahan_c 
      kahan_t = moment0+kahan_y
      kahan_c = (kahan_t-moment0)-kahan_y
      moment0 = kahan_t
!         moment0 = moment0+slip(k, segment)*shear(k, segment)    ! we may need to increase numerical precision
   end do

   call define_slip_field(slip, rake)
   call time_laplace(rupt_time, time_reg)
!
!  An iteration of the simulated annealing algorithm, over each subfault.
!
   subfault = 0
!       begin to perturb       
!
   do k = 1, subfaults
      n_subfault(k) = k
   end do

   do k = 1, subfaults-1
      nran = k
      do while (nran .eq. k .or. nran .gt. subfaults)
         x = ran1()
         nran = int(x*(subfaults-k)+k+1)
      end do
      nbb = n_subfault(nran)
      nn = n_subfault(k)
      n_subfault(k) = nbb
      n_subfault(nran) = nn
   end do

   do k = 1, subfaults
      subfault = n_subfault(k)
      if (subfault .gt. subfaults) stop
      n_total = 0
      do i = 1, segments
         n_total = subfaults_segment(i)+n_total
         if (subfault .le. n_total) then
            segment = i
            subfault_seg = subfault
            exit
         end if
      end do
      do i = 1, segment-1
         subfault_seg = subfault_seg-subfaults_segment(i)
      end do
      d_sub = slip(subfault)
      a_sub = rake(subfault)
      isl = int((t_rise(subfault)-ta0)/dta+0.5)+1
      isr = int((t_fall(subfault)-ta0)/dta+0.5)+1
      rake2 = a_sub*dpi
      a = sin(rake2)*d_sub
      b = cos(rake2)*d_sub
!
!  make up unchange graph
!
      do channel = 1, channels
         delta_freq = delta_freq0/dt_channel(channel)
         omega = -twopi*delta_freq*rupt_time(subfault)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, max_freq
            green_subf = &
            &  (a*green_dip(i, channel, subfault)+b*green_stk(i, channel, subfault)) &
            &  *source(i, channel, isl, isr)*z
            forward_real2(i, channel) = forward_real2(i, channel)-real(green_subf)
            forward_imag2(i, channel) = forward_imag2(i, channel)-aimag(green_subf)
            forward_real(i, channel) = real(forward_real2(i, channel))
            forward_imag(i, channel) = real(forward_imag2(i, channel))
            z = z*z1    ! we may need to increase numuerical precision
         end do
      end do
      if (static) call static_remove_subfault(d_sub, a_sub, subfault)
      if (insar) call insar_remove_subfault(d_sub, a_sub, subfault)
      kahan_y = -slip(subfault)*shear(subfault)-kahan_c 
      kahan_t = moment0+kahan_y
      kahan_c = (kahan_t-moment0)-kahan_y
      moment0 = kahan_t
!  
      n_accept = 0
      npb = np(4*(subfault-1)+1)
      if (npb .lt. 2) exit
!
!  slip extreme values
!
      npb = np(4*(subfault-1)+1)
      dpb = dp(4*(subfault-1)+1)
      slip_beg = beg(4*(subfault-1)+1)
      slip_max = (npb-1)*dpb
      slip_end = slip_beg+slip_max
!
!  rake extreme values
!  
      npb = np(4*(subfault-1)+2)
      dpb = dp(4*(subfault-1)+2)
      angle_beg = beg(4*(subfault-1)+2)
      angle_max = (npb-1)*dpb
      angle_end = angle_beg+angle_max
! 
!  rupture time extreme values.
!
      rupt_beg = time_min(subfault)
      rupt_end = time_max(subfault)
      rupt_max = rupt_end-rupt_beg
!
!  rise time parameters extreme values
!  
      rise_time_beg = ta0
      rise_time_end = ta0+(msou-1)*dta
      rise_time_max = (msou-1)*dta
      do i = 1, max_move
!
!       Save values before the perturbation
!
         t_save = rupt_time(subfault)
         d_save = slip(subfault)
         a_save = rake(subfault)
         l_save = t_rise(subfault)
         r_save = t_fall(subfault)
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
         vuse = rupt_beg - 0.1
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
         misfit2 = 0.d0
!$omp parallel & 
!$omp& default(shared) &
!$omp& private(channel, delta_freq, j, omega, z, z1, green_subf, real1, imag1, coeffs_syn, ex)
!$omp do schedule(static) reduction(+:misfit2)
         do channel = 1, channels
            delta_freq = delta_freq0/dt_channel(channel)
            omega = -twopi*delta_freq*vuse
            z1 = cmplx(cos(omega), sin(omega), double)
            z = cmplx(1.d0, 0.d0, double)
            do j = 1, max_freq
               green_subf = &
               &  (a*green_dip(j, channel, subfault)+b*green_stk(j, channel, subfault)) &
               &  *source(j, channel, isl, isr)*z
               real1(j) = real(green_subf)+forward_real(j, channel)
               imag1(j) = aimag(green_subf)+forward_imag(j, channel)
               z = z*z1    ! we may need to increase numerical precision
            end do
            call wavelet_syn(real1, imag1, coeffs_syn)
            call misfit_channel(channel, coeffs_syn, ex)     
            misfit2 = misfit2 + ex    ! we may need to increase numerical precision
         end do
!$omp end do
!$omp end parallel
         moment0 = moment0+duse*shear(subfault)
         moment = moment0*area
         moment_reg = (moment/moment_input) - 1
         if(abs(moment_reg) .ge. 0.10)then
            moment_reg = sqrt(5*abs(moment_reg)+0.5)
         else
            moment_reg = (10*abs(moment_reg))**4
         endif
!         moment_reg = (moment/moment_input)
         amp = 1.0
         if (static) call static_modify_subfault(duse, ause, subfault, gps_misfit)
         if (insar) call insar_modify_subfault(duse, ause, subfault, insar_misfit)
         call modify_slip_field(subfault, duse, ause)
         call slip_laplace(slip_reg)
         rupt_time(subfault) = vuse
         call time_laplace(rupt_time, time_reg)

         value1 = misfit2+moment_reg*coef_moment+amp*slip_reg*coef_slip
         value1 = value1+coef_time*time_reg
         value1 = value1+coef_gps*gps_misfit+coef_insar*insar_misfit
         moment0 = moment0-duse*shear(subfault)
         diff = value1-current_value
!  
!  Now, we update the kinematic model.
!  
         rand = ran1()
         aux = exp(-diff/t)
         if (aux .gt. rand) then
            current_value = value1
            insar_misfit0 = insar_misfit
            slip(subfault) = duse
            rake(subfault) = ause
            rupt_time(subfault) = vuse
            t_rise(subfault) = (isl-1)*dta+ta0
            t_fall(subfault) = (isr-1)*dta+ta0
            n_accept = n_accept+1
         else
            rupt_time(subfault) = t_save
            t_rise(subfault) = l_save
            t_fall(subfault) = r_save
            slip(subfault) = d_save
            rake(subfault) = a_save
         end if
         min_value = min(min_value, value1)
         if (n_accept .gt. accept_max) exit
      end do
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  finish the perturbation of subevent (segment, subfault_seg)
!
      rake2 = rake(subfault)*dpi
      a = sin(rake2)*slip(subfault)
      b = cos(rake2)*slip(subfault)
      isl = int((t_rise(subfault)-ta0)/dta+0.5)+1
      isr = int((t_fall(subfault)-ta0)/dta+0.5)+1
      kahan_y = slip(subfault)*shear(subfault)-kahan_c 
      kahan_t = moment0+kahan_y
      kahan_c = (kahan_t-moment0)-kahan_y
      moment0 = kahan_t
      d_sub = slip(subfault)
      a_sub = rake(subfault)
      call modify_slip_field(subfault, d_sub, a_sub)
      do channel = 1, channels
         delta_freq = delta_freq0/dt_channel(channel)
         omega = -twopi*delta_freq*rupt_time(subfault)
         z1 = cmplx(cos(omega), sin(omega), double)   
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, max_freq
            green_subf = (a*green_dip(i, channel, subfault)+b*green_stk(i, channel, subfault)) &
            &   *source(i, channel, isl, isr)*z
            forward_real2(i, channel) = forward_real2(i, channel)+real(green_subf)
            forward_imag2(i, channel) = forward_imag2(i, channel)+aimag(green_subf)
            forward_real(i, channel) = real(forward_real2(i, channel))
            forward_imag(i, channel) = real(forward_imag2(i, channel))
            z = z*z1
         end do
      end do
      if (static) call static_add_subfault(d_sub, a_sub, subfault)
      if (insar) call insar_add_subfault(d_sub, a_sub, subfault)
   end do

!
! For insar, now we add ramp perturbation
!
   if (insar) then
      if (present(ramp)) then 
         call insar_remove_ramp(ramp)
         n_accept = 0
!
!  ramp parameters extreme values
!  
         ramp_beg = -50
         ramp_end = 50
         ramp_max = 100
         do i = 1, max_move
!
!       Save values before the perturbation
!
            ramp0(:) = ramp(:)
            ramp1(:) = ramp(:)
!
!  Perturb ramp parameters
!  
            do j = 1, ramp_length
               ramp_use = ramp_beg - 1.
               do while ((ramp_use .lt. ramp_beg) .or. (ramp_use .gt. ramp_end))
                  call cauchy(t, c)
                  ramp_use = ramp1(j)+c*ramp_max
               end do
               ramp1(j) = ramp_use
            end do
            
            call insar_modify_ramp(ramp1, insar_misfit)
            diff = insar_misfit-insar_misfit0
!  
!  Now, we update the ramp.
!  
            rand = ran1()
            aux = exp(-diff/t)
            if (aux .gt. rand) then
               current_value = current_value + coef_insar*(insar_misfit - insar_misfit0)
               insar_misfit0 = insar_misfit
               ramp(:) = ramp1(:)
               n_accept = n_accept+1
            else
               ramp(:) = ramp0(:)
            end if
            if (n_accept .gt. accept_max) exit
         end do
         min_value = min(min_value, current_value)
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  finish the perturbation of ramp
!
         call insar_add_ramp(ramp)
      end if
   end if

   write(*,*) min_value
   end subroutine annealing_iter4


   subroutine annealing_iter5(slip, rake, rupt_time, t_rise, t_fall, t)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  t_rise: array with model risetime values for all subfaults
!  t_fall: array with model falltime values for all subfaults
!  t: Current temperature of the annealing method
!
   use modelling_inputs, only : events, moment_event
   implicit none
   real, intent(inout) :: slip(:), rake(:), rupt_time(:), t_fall(:), t_rise(:)
   real, intent(in) :: t
   integer isl, isr, n_subfault(max_subfaults), n_accept, &
   & nbb, i, k, npb, nn, nran, subfault_seg, segment, channel, subfault, &
   & n_total, j, event
!   real, allocatable :: forward_real(:, :), forward_imag(:, :)
   real :: duse, ause, vuse, &
   & diff, rand, c, aux, dpb, amp, moment_reg(10), moment_reg2, value1, gps_misfit, &
   & moment, d_sub, a_sub, slip_reg, a, b, moment0(10), &
   & time_reg, t_save, a_save, d_save, x, kahan_y, kahan_t, kahan_c, &
   & l_save, r_save, real1(wave_pts2), imag1(wave_pts2), coeffs_syn(wave_pts2), &
   & slip_beg, slip_max, slip_end, angle_beg, angle_end, angle_max, &
   & rupt_beg, rupt_end, rupt_max, rise_time_beg, rise_time_end, rise_time_max
!   real*8, allocatable :: forward_real2(:, :), forward_imag2(:, :)
   real :: delta_freq0, delta_freq, rake2!, ex!, misfit2
   real*8 :: omega, misfit2, ex
   complex :: green_subf
   complex*16 :: z, z1, forward(wave_pts), z0

   z0 = cmplx(0.d0, 0.d0, double)
   value1 = 0.0
   gps_misfit = 0.0
!
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Here, we compute the value of the objective function, 
!  using the input kinematic model.
!
   subfault = 0
   delta_freq0 = 1.0/(2.0**lnpt)

   do channel = 1, channels
      delta_freq = delta_freq0/dt_channel(channel)
      forward(:) = z0
      do subfault = 1, subfaults
         rake2 = rake(subfault)*dpi
         a = sin(rake2)*slip(subfault)
         b = cos(rake2)*slip(subfault)
         isl = int((t_rise(subfault)-ta0)/dta+0.5)+1
         isr = int((t_fall(subfault)-ta0)/dta+0.5)+1
         omega = -twopi*delta_freq*rupt_time(subfault)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, max_freq
!            omega = -twopi*delta_freq*rupt_time(subfault)*(i-1)
!            z = cmplx(cos(omega), sin(omega))
            forward(i) = forward(i) &
            & +(a*green_dip(i, channel, subfault)+b*green_stk(i, channel, subfault)) &
            & *source(i, channel, isl, isr)*z
            z = z*z1    ! we may need to increase numerical precision
         end do
      end do

      do i = 1, wave_pts
         forward_real2(i, channel) = dble(forward(i))
         forward_imag2(i, channel) = aimag(forward(i))
      end do
   end do

   moment_reg2 = 0.0
   do k=1, events
      moment0(k) = 0.0
      kahan_y = 0.0
      kahan_t = 0.0
      kahan_c = 0.0
      do subfault = 1, subfaults
         if (subfault_in_event(subfault, k)) then
            kahan_y = slip(subfault)*shear(subfault)-kahan_c 
            kahan_t = moment0(k)+kahan_y
            kahan_c = (kahan_t-moment0(k))-kahan_y
            moment0(k) = kahan_t
!         moment0 = moment0+slip(k, segment)*shear(k, segment)    ! we may need to increase numerical precision
         endif
      enddo
      moment = moment0(k)*area
      moment_reg(k) = (moment/moment_event(k))-1
      if(abs(moment_reg(k)) .ge. 0.10)then
         moment_reg(k) = sqrt(5*abs(moment_reg(k))+0.5)
      else
         moment_reg(k) = (10*abs(moment_reg(k)))**4
      endif
      moment_reg2 = moment_reg2 + moment_reg(k)
   end do

   call define_slip_field(slip, rake)
   call time_laplace(rupt_time, time_reg)
!
!  An iteration of the simulated annealing algorithm, over each subfault.
!
!       begin to perturb       
!
   do k = 1, subfaults
      n_subfault(k) = k
   end do

   do k = 1, subfaults-1
      nran = k
      do while (nran .eq. k .or. nran .gt. subfaults)
         x = ran1()
         nran = int(x*(subfaults-k)+k+1)
      end do
      nbb = n_subfault(nran)
      nn = n_subfault(k)
      n_subfault(k) = nbb
      n_subfault(nran) = nn
   end do

   do k = 1, subfaults
      subfault = n_subfault(k)
      if (subfault .gt. subfaults) stop
      n_total = 0
      do i = 1, segments
         n_total = subfaults_segment(i)+n_total
         if (subfault .le. n_total) then
            segment = i
            subfault_seg = subfault
            exit
         end if
      end do
      do event=1, events
         if (subfault_in_event(subfault, event)) exit
      enddo
      do i = 1, segment-1
         subfault_seg = subfault_seg-subfaults_segment(i)
      end do
      isl = int((t_rise(subfault)-ta0)/dta+0.5)+1
      isr = int((t_fall(subfault)-ta0)/dta+0.5)+1
      rake2 = rake(subfault)*dpi
      a = sin(rake2)*slip(subfault)
      b = cos(rake2)*slip(subfault)
!
!  make up unchange graph
!
      do channel = 1, channels
         delta_freq = delta_freq0/dt_channel(channel)
         omega = -twopi*delta_freq*rupt_time(subfault)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, max_freq
!            omega = -twopi*delta_freq*(i-1)*rupt_time(subfault)
!            z = cmplx(cos(omega), sin(omega))
            green_subf = &
            & (a*green_dip(i, channel, subfault)+b*green_stk(i, channel, subfault))* &
            & source(i, channel, isl, isr)*z
!            forward_real(i, channel) = forward_real(i, channel)-real(green_subf)
!            forward_imag(i, channel) = forward_imag(i, channel)-aimag(green_subf)
            forward_real2(i, channel) = forward_real2(i, channel)-real(green_subf)
            forward_imag2(i, channel) = forward_imag2(i, channel)-aimag(green_subf)
            forward_real(i, channel) = real(forward_real2(i, channel))
            forward_imag(i, channel) = real(forward_imag2(i, channel))
            z = z*z1    ! we may need to increase numerical precision
         end do
      end do
!      moment0 = moment0-slip(subfault_seg, segment)*shear(subfault_seg, segment)
      kahan_y = -slip(subfault)*shear(subfault)-kahan_c 
      kahan_t = moment0(event)+kahan_y
      kahan_c = (kahan_t-moment0(event))-kahan_y
      moment0(event) = kahan_t
      moment_reg2 = moment_reg2 - moment_reg(event)
!  
      n_accept = 0
      npb = np(4*(subfault-1)+1)
      if (npb .lt. 2) exit
!
!  slip extreme values
!
      npb = np(4*(subfault-1)+1)
      dpb = dp(4*(subfault-1)+1)
      slip_beg = beg(4*(subfault-1)+1)
      slip_max = (npb-1)*dpb
      slip_end = slip_beg+slip_max
!
!  rake extreme values
!  
      npb = np(4*(subfault-1)+2)
      dpb = dp(4*(subfault-1)+2)
      angle_beg = beg(4*(subfault-1)+2)
      angle_max = (npb-1)*dpb
      angle_end = angle_beg+angle_max
! 
!  rupture time extreme values.
!
      rupt_beg = time_min(subfault)
      rupt_end = time_max(subfault)
      rupt_max = rupt_end-rupt_beg
!
!  rise time parameters extreme values
!  
      rise_time_beg = ta0
      rise_time_end = ta0+(msou-1)*dta
      rise_time_max = (msou-1)*dta
      do i = 1, max_move
!
!       Save values before the perturbation
!
         t_save = rupt_time(subfault)
         d_save = slip(subfault)
         a_save = rake(subfault)
         l_save = t_rise(subfault)
         r_save = t_fall(subfault)
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
         vuse = rupt_beg - 0.1
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
         misfit2 = 0.d0
!$omp parallel & 
!$omp& default(shared) &
!$omp& private(channel, delta_freq, j, omega, z, z1, green_subf, real1, imag1, coeffs_syn, ex)
!$omp do schedule(static) reduction(+:misfit2)
         do channel = 1, channels
            delta_freq = delta_freq0/dt_channel(channel)
            omega = -twopi*delta_freq*vuse
            z1 = cmplx(cos(omega), sin(omega), double) 
            z = cmplx(1.d0, 0.d0, double)
            do j = 1, max_freq
!               omega = -twopi_0*delta_freq*(i-1)*vuse
!               z = cmplx(cos(omega), sin(omega)) 
               green_subf = &
               & (a*green_dip(j, channel, subfault)+b*green_stk(j, channel, subfault))* &
               & source(j, channel, isl, isr)*z
               real1(j) = real(green_subf)+forward_real(j, channel)
               imag1(j) = aimag(green_subf)+forward_imag(j, channel)
               z = z*z1    ! we may need to increase numerical precision
            end do
            call wavelet_syn(real1, imag1, coeffs_syn)
            call misfit_channel(channel, coeffs_syn, ex)     
            misfit2 = misfit2 + ex    ! we may need to increase numerical precision
         end do
!$omp end do
!$omp end parallel
         moment0(event) = moment0(event)+duse*shear(subfault)
         moment = moment0(event)*area
         moment_reg(event) = (moment/moment_event(event))-1
         if(abs(moment_reg(event)) .ge. 0.10)then
            moment_reg(event) = sqrt(5*abs(moment_reg(event))+0.5)
         else
            moment_reg(event) = (10*abs(moment_reg(event)))**4
         endif
         moment_reg2 = moment_reg(event) + moment_reg2
!         moment_reg = (moment/moment_input)
         amp = 1.0
         call modify_slip_field(subfault, duse, ause)
         call slip_laplace(slip_reg)
         rupt_time(subfault) = vuse
         call time_laplace(rupt_time, time_reg)

         value1 = misfit2+moment_reg2*coef_moment+amp*slip_reg*coef_slip
         value1 = value1+coef_time*time_reg+coef_gps*gps_misfit
         moment0(event) = moment0(event)-duse*shear(subfault)
         moment_reg2 = moment_reg2 - moment_reg(event)
         diff = value1-current_value
!  
!  Now, we update the kinematic model.
!  
         rand = ran1()
         aux = exp(-diff/t)
         if (aux .gt. rand) then
            current_value = value1
            slip(subfault) = duse
            rake(subfault) = ause
            rupt_time(subfault) = vuse
            t_rise(subfault) = (isl-1)*dta+ta0
            t_fall(subfault) = (isr-1)*dta+ta0
            n_accept = n_accept+1
         else
            rupt_time(subfault) = t_save
            t_rise(subfault) = l_save
            t_fall(subfault) = r_save
            slip(subfault) = d_save
            rake(subfault) = a_save
         end if
         min_value = min(min_value, value1)
         if (n_accept .gt. accept_max) exit
      end do
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  finish the perturbation of subevent (segment, subfault_seg)
!
      rake2 = rake(subfault)*dpi
      a = sin(rake2)*slip(subfault)
      b = cos(rake2)*slip(subfault)
      isl = int((t_rise(subfault)-ta0)/dta+0.5)+1
      isr = int((t_fall(subfault)-ta0)/dta+0.5)+1
      kahan_y = slip(subfault)*shear(subfault)-kahan_c 
      kahan_t = moment0(event)+kahan_y
      kahan_c = (kahan_t-moment0(event))-kahan_y
      moment0(event) = kahan_t
      moment = moment0(event)*area
      moment_reg(event) = (moment/moment_event(event))-1
      if(abs(moment_reg(event)) .ge. 0.10)then
         moment_reg(event) = sqrt(5*abs(moment_reg(event))+0.5)
      else
         moment_reg(event) = (10*abs(moment_reg(event)))**4
      endif
      moment_reg2 = moment_reg2 + moment_reg(event)
      d_sub = slip(subfault)
      a_sub = rake(subfault)
      call modify_slip_field(subfault, d_sub, a_sub)
      do channel = 1, channels
         delta_freq = delta_freq0/dt_channel(channel)
         omega = -twopi*delta_freq*rupt_time(subfault)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, max_freq
            green_subf = (a*green_dip(i, channel, subfault)+b*green_stk(i, channel, subfault)) &
            &   *source(i, channel, isl, isr)*z
            forward_real2(i, channel) = forward_real2(i, channel)+real(green_subf)
            forward_imag2(i, channel) = forward_imag2(i, channel)+aimag(green_subf)
            forward_real(i, channel) = real(forward_real2(i, channel))
            forward_imag(i, channel) = real(forward_imag2(i, channel))
            z = z*z1
         end do
      end do

   end do
   write(*,*) min_value
   end subroutine annealing_iter5 

   
   subroutine annealing_iter6(slip, rake, rupt_time, t_rise, &
   & t_fall, t, static, insar, ramp)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  t_rise: array with model risetime values for all subfaults
!  t_fall: array with model falltime values for all subfaults
!  t: Current temperature of the annealing method
!  static: True if static GPS data used in modelling, False otherwise
!  insar: True if insar data used in modelling, False otherwise
!  ramp: Value of insar ramp, optional
!
   use modelling_inputs, only : events, moment_event
   implicit none
   real, intent(inout) :: slip(:), rake(:), rupt_time(:), t_fall(:), t_rise(:)
   real, intent(in) :: t
   logical, intent(in) :: static, insar
   real*8, optional :: ramp(:)
   integer isl, isr, n_subfault(max_subfaults), n_accept, &
   & nbb, i, k, npb, nn, nran, subfault_seg, segment, channel, subfault, &
   & n_total, j, event
!   real, allocatable :: forward_real(:, :), forward_imag(:, :)
   real :: duse, ause, vuse, &
   & diff, rand, c, aux, dpb, amp, moment_reg(10), moment_reg2, value1, gps_misfit, insar_misfit, &
   & moment, d_sub, a_sub, slip_reg, a, b, kahan_y, kahan_c, kahan_t, &
   & time_reg, t_save, a_save, d_save, x, moment0(10), &
   & l_save, r_save, real1(wave_pts2), imag1(wave_pts2), coeffs_syn(wave_pts2), &
   & slip_beg, slip_max, slip_end, angle_beg, angle_end, angle_max, &
   & rupt_beg, rupt_end, rupt_max, rise_time_beg, rise_time_end, rise_time_max
   real ramp_beg, ramp_end, ramp_max, ramp_use
   real*8 :: ramp0(36), ramp1(36)
!   real*8, allocatable :: forward_real2(:, :), forward_imag2(:, :)
   real*8 :: omega, misfit2, ex
   real :: delta_freq, delta_freq0, rake2!, ex
   complex :: green_subf
   complex*16 :: z, z1, forward(wave_pts), z0
!
   z0 = cmplx(0.d0, 0.d0, double)
   value1 = 0.0
!
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Here, we compute the value of the objective function, 
!  using the input kinematic model.
!
   subfault = 0
   delta_freq0 = 1.0/(2.0**lnpt)

   do channel = 1, channels
      delta_freq = delta_freq0/dt_channel(channel)
      forward(:) = z0
      do subfault = 1, subfaults
         rake2 = rake(subfault)*dpi
         a = sin(rake2)*slip(subfault)
         b = cos(rake2)*slip(subfault)
         isl = int((t_rise(subfault)-ta0)/dta+0.5)+1
         isr = int((t_fall(subfault)-ta0)/dta+0.5)+1
         omega = -twopi*delta_freq*rupt_time(subfault)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, max_freq
            forward(i) = forward(i) &
            & +(a*green_dip(i, channel, subfault)+b*green_stk(i, channel, subfault)) &
            & *source(i, channel, isl, isr)*z
            z = z*z1    ! we may need to increase numerical precision
         end do
      end do

      do i = 1, wave_pts
         forward_real2(i, channel) = dble(forward(i))
         forward_imag2(i, channel) = aimag(forward(i))
      end do
   end do
   if (static) call static_synthetic(slip, rake, gps_misfit)
   if (insar) then
      if (present(ramp)) then
         call insar_synthetic(slip, rake, insar_misfit, ramp)
      else
         call insar_synthetic(slip, rake, insar_misfit)
      endif
   endif

   moment_reg2 = 0.0
   do k=1, events
      moment0(k) = 0.0
      kahan_y = 0.0
      kahan_t = 0.0
      kahan_c = 0.0
      do subfault = 1, subfaults
         if (subfault_in_event(subfault, k)) then
            kahan_y = slip(subfault)*shear(subfault)-kahan_c 
            kahan_t = moment0(k)+kahan_y
            kahan_c = (kahan_t-moment0(k))-kahan_y
            moment0(k) = kahan_t
!         moment0 = moment0+slip(k, segment)*shear(k, segment)    ! we may need to increase numerical precision
         endif
      enddo
      moment = moment0(k)*area
      moment_reg(k) = (moment/moment_event(k))-1
      if(abs(moment_reg(k)) .ge. 0.10)then
         moment_reg(k) = sqrt(5*abs(moment_reg(k))+0.5)
      else
         moment_reg(k) = (10*abs(moment_reg(k)))**4
      endif
      moment_reg2 = moment_reg2 + moment_reg(k)
   end do

   call define_slip_field(slip, rake)
   call time_laplace(rupt_time, time_reg)
!
!  An iteration of the simulated annealing algorithm, over each subfault.
!
   subfault = 0
!       begin to perturb       
!
   do k = 1, subfaults
      n_subfault(k) = k
   end do

   do k = 1, subfaults-1
      nran = k
      do while (nran .eq. k .or. nran .gt. subfaults)
         x = ran1()
         nran = int(x*(subfaults-k)+k+1)
      end do
      nbb = n_subfault(nran)
      nn = n_subfault(k)
      n_subfault(k) = nbb
      n_subfault(nran) = nn
   end do

   do k = 1, subfaults
      subfault = n_subfault(k)
      if (subfault .gt. subfaults) stop
      n_total = 0
      do i = 1, segments
         n_total = subfaults_segment(i)+n_total
         if (subfault .le. n_total) then
            segment = i
            subfault_seg = subfault
            exit
         end if
      end do
      do event=1, events
         if (subfault_in_event(subfault, event)) exit
      enddo
      do i = 1, segment-1
         subfault_seg = subfault_seg-subfaults_segment(i)
      end do
      d_sub = slip(subfault)
      a_sub = rake(subfault)
      isl = int((t_rise(subfault)-ta0)/dta+0.5)+1
      isr = int((t_fall(subfault)-ta0)/dta+0.5)+1
      rake2 = a_sub*dpi
      a = sin(rake2)*d_sub
      b = cos(rake2)*d_sub
!
!  make up unchange graph
!
      do channel = 1, channels
         delta_freq = delta_freq0/dt_channel(channel)
         omega = -twopi*delta_freq*rupt_time(subfault)
         z1 = cmplx(cos(omega), sin(omega), double)
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, max_freq
            green_subf = &
            &  (a*green_dip(i, channel, subfault)+b*green_stk(i, channel, subfault)) &
            &  *source(i, channel, isl, isr)*z
            forward_real2(i, channel) = forward_real2(i, channel)-real(green_subf)
            forward_imag2(i, channel) = forward_imag2(i, channel)-aimag(green_subf)
            forward_real(i, channel) = real(forward_real2(i, channel))
            forward_imag(i, channel) = real(forward_imag2(i, channel))
            z = z*z1    ! we may need to increase numuerical precision
         end do
      end do
      if (static) call static_remove_subfault(d_sub, a_sub, subfault)
      if (insar) call insar_remove_subfault(d_sub, a_sub, subfault)
      kahan_y = -slip(subfault)*shear(subfault)-kahan_c 
      kahan_t = moment0(event)+kahan_y
      kahan_c = (kahan_t-moment0(event))-kahan_y
      moment0(event) = kahan_t
      moment_reg2 = moment_reg2 - moment_reg(event)
!  
      n_accept = 0
      npb = np(4*(subfault-1)+1)
      if (npb .lt. 2) exit
!
!  slip extreme values
!
      npb = np(4*(subfault-1)+1)
      dpb = dp(4*(subfault-1)+1)
      slip_beg = beg(4*(subfault-1)+1)
      slip_max = (npb-1)*dpb
      slip_end = slip_beg+slip_max
!
!  rake extreme values
!  
      npb = np(4*(subfault-1)+2)
      dpb = dp(4*(subfault-1)+2)
      angle_beg = beg(4*(subfault-1)+2)
      angle_max = (npb-1)*dpb
      angle_end = angle_beg+angle_max
! 
!  rupture time extreme values.
!
      rupt_beg = time_min(subfault)
      rupt_end = time_max(subfault)
      rupt_max = rupt_end-rupt_beg
!
!  rise time parameters extreme values
!  
      rise_time_beg = ta0
      rise_time_end = ta0+(msou-1)*dta
      rise_time_max = (msou-1)*dta
      do i = 1, max_move
!
!       Save values before the perturbation
!
         t_save = rupt_time(subfault)
         d_save = slip(subfault)
         a_save = rake(subfault)
         l_save = t_rise(subfault)
         r_save = t_fall(subfault)
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
         vuse = rupt_beg - 0.1
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
         misfit2 = 0.d0
!$omp parallel & 
!$omp& default(shared) &
!$omp& private(channel, delta_freq, j, omega, z, z1, green_subf, real1, imag1, coeffs_syn, ex)
!$omp do schedule(static) reduction(+:misfit2)
         do channel = 1, channels
            delta_freq = delta_freq0/dt_channel(channel)
            omega = -twopi*delta_freq*vuse
            z1 = cmplx(cos(omega), sin(omega), double)
            z = cmplx(1.d0, 0.d0, double)
            do j = 1, max_freq
               green_subf = &
               &  (a*green_dip(j, channel, subfault)+b*green_stk(j, channel, subfault)) &
               &  *source(j, channel, isl, isr)*z
               real1(j) = real(green_subf)+forward_real(j, channel)
               imag1(j) = aimag(green_subf)+forward_imag(j, channel)
               z = z*z1    ! we may need to increase numerical precision
            end do
            call wavelet_syn(real1, imag1, coeffs_syn)
            call misfit_channel(channel, coeffs_syn, ex)     
            misfit2 = misfit2 + ex    ! we may need to increase numerical precision
         end do
!$omp end do
!$omp end parallel
         moment0(event) = moment0(event)+duse*shear(subfault)
         moment = moment0(event)*area
         moment_reg(event) = (moment/moment_event(event))-1
         if(abs(moment_reg(event)) .ge. 0.10)then
            moment_reg(event) = sqrt(5*abs(moment_reg(event))+0.5)
         else
            moment_reg(event) = (10*abs(moment_reg(event)))**4
         endif
         moment_reg2 = moment_reg(event) + moment_reg2
!         moment_reg = (moment/moment_input)
         amp = 1.0
         if (static) call static_modify_subfault(duse, ause, subfault, gps_misfit)
         if (insar) call insar_modify_subfault(duse, ause, subfault, insar_misfit)
         call modify_slip_field(subfault, duse, ause)
         call slip_laplace(slip_reg)
         rupt_time(subfault) = vuse
         call time_laplace(rupt_time, time_reg)

         value1 = misfit2+moment_reg2*coef_moment+amp*slip_reg*coef_slip
         value1 = value1+coef_time*time_reg
         value1 = value1+coef_gps*gps_misfit+coef_insar*insar_misfit
         moment0(event) = moment0(event)-duse*shear(subfault)
         moment_reg2 = moment_reg2 - moment_reg(event)
         diff = value1-current_value
!  
!  Now, we update the kinematic model.
!  
         rand = ran1()
         aux = exp(-diff/t)
         if (aux .gt. rand) then
            current_value = value1
            insar_misfit0 = insar_misfit
            slip(subfault) = duse
            rake(subfault) = ause
            rupt_time(subfault) = vuse
            t_rise(subfault) = (isl-1)*dta+ta0
            t_fall(subfault) = (isr-1)*dta+ta0
            n_accept = n_accept+1
         else
            rupt_time(subfault) = t_save
            t_rise(subfault) = l_save
            t_fall(subfault) = r_save
            slip(subfault) = d_save
            rake(subfault) = a_save
         end if
         min_value = min(min_value, value1)
         if (n_accept .gt. accept_max) exit
      end do
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  finish the perturbation of subevent (segment, subfault_seg)
!
      rake2 = rake(subfault)*dpi
      a = sin(rake2)*slip(subfault)
      b = cos(rake2)*slip(subfault)
      isl = int((t_rise(subfault)-ta0)/dta+0.5)+1
      isr = int((t_fall(subfault)-ta0)/dta+0.5)+1
      kahan_y = slip(subfault)*shear(subfault)-kahan_c 
      kahan_t = moment0(event)+kahan_y
      kahan_c = (kahan_t-moment0(event))-kahan_y
      moment0(event) = kahan_t
      moment = moment0(event)*area
      moment_reg(event) = (moment/moment_event(event))-1
      if(abs(moment_reg(event)) .ge. 0.10)then
         moment_reg(event) = sqrt(5*abs(moment_reg(event))+0.5)
      else
         moment_reg(event) = (10*abs(moment_reg(event)))**4
      endif
      moment_reg2 = moment_reg(event) + moment_reg2
      d_sub = slip(subfault)
      a_sub = rake(subfault)
      call modify_slip_field(subfault, d_sub, a_sub)
      do channel = 1, channels
         delta_freq = delta_freq0/dt_channel(channel)
         omega = -twopi*delta_freq*rupt_time(subfault)
         z1 = cmplx(cos(omega), sin(omega), double)   
         z = cmplx(1.d0, 0.d0, double)
         do i = 1, max_freq
            green_subf = (a*green_dip(i, channel, subfault)+b*green_stk(i, channel, subfault)) &
            &   *source(i, channel, isl, isr)*z
            forward_real2(i, channel) = forward_real2(i, channel)+real(green_subf)
            forward_imag2(i, channel) = forward_imag2(i, channel)+aimag(green_subf)
            forward_real(i, channel) = real(forward_real2(i, channel))
            forward_imag(i, channel) = real(forward_imag2(i, channel))
            z = z*z1
         end do
      end do
      if (static) call static_add_subfault(d_sub, a_sub, subfault)
      if (insar) call insar_add_subfault(d_sub, a_sub, subfault)
   end do

!
! For insar, now we add ramp perturbation
!
   if (insar) then
      if (present(ramp)) then 
         call insar_remove_ramp(ramp)
         n_accept = 0
!
!  ramp parameters extreme values
!  
         ramp_beg = -50
         ramp_end = 50
         ramp_max = 100
         do i = 1, max_move
!
!       Save values before the perturbation
!
            ramp0(:) = ramp(:)
            ramp1(:) = ramp(:)
!
!  Perturb ramp parameters
!  
            do j = 1, ramp_length
               ramp_use = ramp_beg - 1.
               do while ((ramp_use .lt. ramp_beg) .or. (ramp_use .gt. ramp_end))
                  call cauchy(t, c)
                  ramp_use = ramp1(j)+c*ramp_max
               end do
               ramp1(j) = ramp_use
            end do
            
            call insar_modify_ramp(ramp1, insar_misfit)
            diff = insar_misfit-insar_misfit0
!  
!  Now, we update the ramp.
!  
            rand = ran1()
            aux = exp(-diff/t)
            if (aux .gt. rand) then
               current_value = current_value + coef_insar*(insar_misfit - insar_misfit0)
               insar_misfit0 = insar_misfit
               ramp(:) = ramp1(:)
               n_accept = n_accept+1
            else
               ramp(:) = ramp0(:)
            end if
            if (n_accept .gt. accept_max) exit
         end do
         min_value = min(min_value, current_value)
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  finish the perturbation of ramp
!
         call insar_add_ramp(ramp)
      end if
   end if

   write(*,*) min_value
   end subroutine annealing_iter6

end module annealing 
