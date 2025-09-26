module annealing_static


   use constants, only : pi, max_stations, dpi, twopi, max_subf, max_seg, wave_pts2, &
            &   max_subfaults2, wave_pts, max_subfaults
   use random_gen, only : ran1, cauchy
   use modelling_inputs, only : get_weights_moment_end 
   use regularization, only : slip_laplace, define_slip_field, modify_slip_field
   use static_data, only : static_synthetic, static_remove_subfault, &
                       &   static_modify_subfault, static_add_subfault
   use imagery_data, only : imagery_synthetic, imagery_remove_subfault, &
                       &   imagery_modify_subfault, imagery_add_subfault, &
                       &   imagery_remove_ramp, imagery_modify_ramp, &
                       &   imagery_add_ramp, ramp_length
   use omp_lib
   implicit none
   real :: coef_moment, coef_slip, coef_gnss, coef_imagery
   real :: current_value, min_value, min_dt, area
   real :: imagery_misfit0
   real :: moment_input, smooth_moment, smooth_slip, smooth_time, emin0
   integer :: subfaults_segment(max_seg)
   integer, parameter :: double = kind(1.d0)
   integer, private :: threads
   integer, parameter, private :: max_move=50, accept_max=5
   real :: minimum(max_subfaults2), delta(max_subfaults2)
   integer :: n_values(max_subfaults2), segments, subfaults
   real :: shear(max_subfaults), dxs, dys
   integer :: nxs_sub(max_seg), nys_sub(max_seg) 


contains


   subroutine annealingstatic_set_fault_properties()
   use model_parameters, only : get_shear, get_segments, get_subfaults, get_space
   implicit none
   real :: dip(max_seg), strike(max_seg), delay_seg(max_seg)
   integer :: cum_subfaults(max_seg), nx_p, ny_p
   real :: v_min, v_max, v_ref, time_ref(max_subfaults2)
   real :: time_min(max_subfaults2), time_max(max_subfaults2)
   call get_shear(shear)
   call get_segments(nxs_sub, nys_sub, dip, strike, delay_seg, segments, subfaults, cum_subfaults)
   call get_subfaults(dxs, dys, nx_p, ny_p, v_min, v_max, v_ref)
   call get_space(time_min, time_max, time_ref, minimum, delta, n_values)
   end subroutine annealingstatic_set_fault_properties
   

   subroutine annealingstatic_set_procedure_param()
   implicit none
   real :: real0
   call get_weights_moment_end(moment_input, smooth_moment, smooth_slip, smooth_time, real0, emin0)
   end subroutine annealingstatic_set_procedure_param


   subroutine print_static_summary(slip, rake, static, imagery, get_coeff, ramp)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  static: True if static GNSS data used in modelling, False otherwise
!  imagery: True if imagery data used in modelling, False otherwise
!  get_coeff: get regularization coefficients if and only if this is True
!  ramp: Value of imagery ramp, optional
!
   implicit none
   real*8, optional :: ramp(:)
   real :: slip(:), rake(:)
   real amp, moment, moment_reg, dt, value1, er0, slip_reg, gnss_misfit, imagery_misfit
   real :: delta_freq, delta_freq0, moment0, kahan_y, kahan_t, kahan_c
   integer :: i, segment, channel, isl, isr, ixs, iys, jf, k, subfault, subfault_seg 
   real*8 :: misfit2
   logical :: static, get_coeff, imagery

   do segment = 1, segments
      subfaults_segment(segment) = nys_sub(segment)*nxs_sub(segment)
   end do
   gnss_misfit = 0.0
   imagery_misfit = 0.0
!
! Compute synthetics given current fault model
!
   area = dxs*dys*(1.e+10)

   amp = 1.0
   moment0 = 0.0
   kahan_y = 0.0
   kahan_t = 0.0
   kahan_c = 0.0
   subfault = 0
   do subfault = 1, subfaults
      kahan_y = slip(subfault)*shear(subfault)-kahan_c 
      kahan_t = moment0+kahan_y
      kahan_c = (kahan_t-moment0)-kahan_y
      moment0 = kahan_t
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
   call define_slip_field(slip, rake)
   call slip_laplace(slip_reg)
 
   coef_gnss = 0.0
   coef_imagery = 0.0
   imagery_misfit0 = 0.0
   if (static) call static_synthetic(slip, rake, gnss_misfit)
   if (imagery) then
      if (present(ramp)) then
         call imagery_synthetic(slip, rake, imagery_misfit, ramp)
      else
         call imagery_synthetic(slip, rake, imagery_misfit)
      endif
      imagery_misfit0 = imagery_misfit
   endif
   misfit2 = imagery_misfit + gnss_misfit

   if (get_coeff) then
      coef_moment=smooth_moment*(misfit2 - current_value)
      coef_moment = min(coef_moment, 1.0)
      coef_slip = smooth_slip*misfit2/(slip_reg*amp)
      coef_slip = min(0.003, coef_slip)
!      coef_time = smooth_time*misfit2/(time_reg*amp)
!      if (static) then       !! update!
!         coef_gnss = misfit2/(gnss_misfit*amp)
!      endif
!      if (imagery) then
!         coef_imagery = misfit2/(imagery_misfit*amp)
!      endif
   endif

   value1 = misfit2 + coef_moment*moment_reg+coef_slip*slip_reg*amp
!   value1 = value1+coef_gnss*gnss_misfit+coef_imagery*imagery_misfit
   write(*,'()')
   write(*,*)'averaged misfit error', misfit2
   write(*,*)'moment error', moment_reg
   write(*,*)'slip smoothness penalization', slip_reg
   if (static) write(*,*)'static data misfit', gnss_misfit
   if (imagery) write(*,*)'imagery data misfit', imagery_misfit
   write(*,*)'total moment of the inversion', moment
   write(*,'(/A, F10.7)')'moment error coefficient', coef_moment
   write(*,*)'slip smoothness penalization coefficient', coef_slip
!   if (static) write(*,*)'static data penalization coefficient', coef_gnss
!   if (imagery) write(*,*)'imagery data penalization coefficient', coef_imagery
   write(*,'(/A, I4)')'Amount of variables: ', 5 * subfaults
!   write(*,*)'Amount of data values: ', used_data
   current_value = value1
   min_value = value1
   open(12,file='modelling_summary.txt')
   write(12,'(/A/)')'Modelling Report'
   write(12,*)'averaged misfit error', misfit2
   write(12,*)'moment error', moment_reg
   write(12,*)'slip smoothness penalization', slip_reg
   if (static) write(12,*)'static data misfit', gnss_misfit
   if (imagery) write(12,*)'imagery data misfit', imagery_misfit
   write(12,*)'objective function value', value1
   write(12,*)'total moment of the inversion', moment
   write(12,'(/A, F10.7)')'moment error coefficient', coef_moment
   write(12,*)'slip smoothness penalization coefficient', coef_slip
!   if (static) write(12,*)'static data penalization coefficient', coef_gnss
!   if (imagery) write(12,*)'imagery data penalization coefficient', coef_imagery
   write(12,'(/A, I4)')'Amount of variables: ', 5 * subfaults
!   write(12,*)'Amount of data values: ', used_data
   close(12)
   end subroutine print_static_summary

   
   subroutine annealing_iter(slip, rake, t, static, imagery, ramp)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  t: Current temperature of the annealing method
!  static: True if static GNSS data used in modelling, False otherwise
!  imagery: True if imagery data used in modelling, False otherwise
!  ramp: Value of imagery ramp, optional
!
   implicit none
   integer isl, isr, n_subfault(max_subfaults), n_accept, &
   & subfault0, i, k, values0, subfault1, index0, subfault_seg, segment, channel, subfault, iys, &
   & ixs, n_total, j
   real*8, optional :: ramp(:)
   real slip(:), rake(:), t, slip1, rake1, &
   & diff, random0, aux, delta0, amp, moment_reg, value1, gnss_misfit, imagery_misfit, &
   & moment, slip0, rake0, slip_reg, kahan_y, kahan_c, kahan_t, &
   & time_reg, moment0, &
   & slip_beg, slip_max, slip_end, angle_beg, angle_end, angle_max
   real ramp_beg, ramp_end, ramp_max, ramp_use
   real*8 :: ramp0(36), ramp1(36)
   real*8 :: omega, misfit2, ex
   real :: delta_freq, delta_freq0!, ex
   logical :: static, imagery
!
   value1 = 0.0
   gnss_misfit = 0.0
   imagery_misfit = 0.0
!
!  ++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Here, we compute the value of the objective function, 
!  using the input kinematic model.
!
   subfault = 0

   if (static) call static_synthetic(slip, rake, gnss_misfit)
   if (imagery) then
      if (present(ramp)) then
         call imagery_synthetic(slip, rake, imagery_misfit, ramp)
      else
         call imagery_synthetic(slip, rake, imagery_misfit)
      endif
   endif
   misfit2 = imagery_misfit + gnss_misfit
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
      index0 = k
      do while (index0 .eq. k .or. index0 .gt. subfaults)
         !call random_number(x)
         random0 = ran1()
         index0 = int(random0*(subfaults-k)+k+1)
      end do
      subfault0 = n_subfault(index0)
      subfault1 = n_subfault(k)
      n_subfault(k) = subfault0
      n_subfault(index0) = subfault1
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
      slip0 = slip(subfault)
      rake0 = rake(subfault)
!
!  make up unchange graph
!
      if (static) call static_remove_subfault(slip0, rake0, subfault)
      if (imagery) call imagery_remove_subfault(slip0, rake0, subfault)
      kahan_y = -slip(subfault)*shear(subfault)-kahan_c 
      kahan_t = moment0+kahan_y
      kahan_c = (kahan_t-moment0)-kahan_y
      moment0 = kahan_t
!  
      n_accept = 0
      values0 = n_values(4*(subfault-1)+1)
      if (values0 .lt. 2) exit
!
!  slip extreme values
!
      values0 = n_values(4*(subfault-1)+1)
      delta0 = delta(4*(subfault-1)+1)
      slip_beg = minimum(4*(subfault-1)+1)
      slip_max = (values0-1)*delta0
      slip_end = slip_beg+slip_max
!
!  rake extreme values
!  
      values0 = n_values(4*(subfault-1)+2)
      delta0 = delta(4*(subfault-1)+2)
      angle_beg = minimum(4*(subfault-1)+2)
      angle_max = (values0-1)*delta0
      angle_end = angle_beg+angle_max
      do i = 1, max_move
!
!       Save values before the perturbation
!
         slip0 = slip(subfault)
         rake0 = rake(subfault)
!
!  Perturb the slip
!
         slip1 = slip_beg - 1.
         do while ((slip1 .le. slip_beg) .or. (slip1 .ge. slip_end))
            call cauchy(t, random0)                           
            slip1 = slip0+random0*slip_max
         end do
!
!  Perturb the rake
!
         rake1 = angle_beg - 1.
         do while ((rake1 .lt. angle_beg) .or. (rake1 .gt. angle_end))
            call cauchy(t, random0)                          
            rake1 = rake0+random0*angle_max
         end do
         
         moment0 = moment0+slip1*shear(subfault)
         moment = moment0*area
         moment_reg = (moment/moment_input) - 1
         if(abs(moment_reg) .ge. 0.10)then
            moment_reg = sqrt(5*abs(moment_reg)+0.5)
         else
            moment_reg = (10*abs(moment_reg))**4
         endif
!         moment_reg = (moment/moment_input)
         amp = 1.0
         if (static) call static_modify_subfault(slip1, rake1, subfault, gnss_misfit)
         if (imagery) call imagery_modify_subfault(slip1, rake1, subfault, imagery_misfit)
         call modify_slip_field(subfault, slip1, rake1)
         call slip_laplace(slip_reg)

         misfit2 = 0.d0
         if (static) misfit2 = misfit2 + gnss_misfit
         if (imagery) misfit2 = misfit2 + imagery_misfit
!         misfit2 = imagery_misfit + gnss_misfit
         value1 = misfit2 + moment_reg*coef_moment+amp*slip_reg*coef_slip
!         value1 = value1+coef_gnss*gnss_misfit+coef_imagery*imagery_misfit
         moment0 = moment0-slip1*shear(subfault)
         diff = value1-current_value
!  
!  Now, we update the kinematic model.
! 
         !call random_number(rand) 
         random0 = ran1()
         aux = exp(-diff/t)
         if (aux .gt. random0) then
            current_value = value1
            imagery_misfit0 = imagery_misfit
            slip(subfault) = slip1
            rake(subfault) = rake1
            n_accept = n_accept+1
         else
            slip(subfault) = slip0
            rake(subfault) = rake0
         end if
         min_value = min(min_value, value1)
         if (n_accept .gt. accept_max) exit
      end do
!
!  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  finish the perturbation of subevent (segment, subfault_seg)
!
      kahan_y = slip(subfault)*shear(subfault)-kahan_c 
      kahan_t = moment0+kahan_y
      kahan_c = (kahan_t-moment0)-kahan_y
      moment0 = kahan_t
      slip0 = slip(subfault)
      rake0 = rake(subfault)
      call modify_slip_field(subfault, slip0, rake0)
      if (static) call static_add_subfault(slip0, rake0, subfault)
      if (imagery) call imagery_add_subfault(slip0, rake0, subfault)
   end do

!
! For imagery, now we add ramp perturbation
!
   if (imagery) then
      if (present(ramp)) then 
         call imagery_remove_ramp(ramp)
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
                  call cauchy(t, random0)
                  ramp_use = ramp1(j)+random0*ramp_max
               end do
               ramp1(j) = ramp_use
            end do
            
            call imagery_modify_ramp(ramp1, imagery_misfit)
            diff = imagery_misfit-imagery_misfit0
!  
!  Now, we update the ramp.
! 
            !call random_number(rand) 
            random0 = ran1()
            aux = exp(-diff/t)
            if (aux .gt. random0) then
               current_value = current_value + (imagery_misfit - imagery_misfit0)
               imagery_misfit0 = imagery_misfit
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
         call imagery_add_ramp(ramp)
      end if
   end if

   write(*,*) min_value
   end subroutine annealing_iter


end module annealing_static 
