module ffm_methods


   use constants, only : max_seg, max_subfaults, max_subf
   use model_parameters, only : write_model, deallocate_ps
   use modelling_inputs, only : get_annealing_param
   use get_stations_data, only : get_data
   use retrieve_gf, only : get_gf, deallocate_gf, retrievegf_set_data_properties, &
                    &   retrievegf_set_fault_parameters
   use wavelets, only : wavelets_set_data_properties, fourier_coefs, meyer_yamada
   use misfit_eval, only : misfit_eval_set_data_properties
   use save_forward, only : write_forward, saveforward_set_fault_parameters, &
                    &   saveforward_set_data_properties
   use rise_time, only : get_source_fun, deallocate_source
   use static_data, only : initial_gnss, staticdata_set_fault_parameters
   use insar_data, only : initial_insar, get_insar_gf, deallocate_insar_gf, &
                    &   get_insar_data, is_ramp, initial_ramp, &
                    &   insardata_set_fault_parameters
   use annealing, only : initial_model, print_summary, print_summary2, &
                    &   annealing_iter3, annealing_iter4, n_threads, &
                    &   annealing_iter5, annealing_iter6, &
                    &   annealing_set_data_properties, annealing_set_fault_parameters, &
                    &   annealing_set_procedure_param, annealing_set_events, &
                    &   allocate_forward, deallocate_forward
   use annealing_static, only : print_static_summary, annealing_iter, &
                    &   annealingstatic_set_fault_properties, &
                    &   annealingstatic_set_procedure_param
   implicit none
   real :: rupt_time0(max_subfaults)
   real :: t_rise0(max_subfaults), t_fall0(max_subfaults)
   real :: t, t0, cooling_rate, t_stop, t_mid
   real*8 :: ramp(36)
   integer :: i, n_iter, start_annealing
   character(len=10) :: input


contains
   

   subroutine ffmmethod_set_procedure_param()
!
!  Args:
!
   implicit none
   integer :: seed0
   call get_annealing_param(n_iter, seed0, t0, cooling_rate, t_stop, start_annealing, t_mid)
   end subroutine ffmmethod_set_procedure_param


   subroutine check_waveforms(strong, cgnss, body, surf, dart, use_waveforms)
!
!  Args:
!  strong: True if strong motion data are used, False otherwise
!  cgnss: True if cGNSS data are used, False otherwise
!  body: True if body wave data are used, False otherwise
!  surf: True if surface wave data are used, False otherwise
!  dart: True if DART data are used, False otherwise
!  use_waveforms: True if some waveform data used in modelling, False otherwise
!
   implicit none
   logical :: strong, cgnss, dart, body, surf, use_waveforms
   use_waveforms = .False.
   if (strong) use_waveforms = .True.
   if (cgnss) use_waveforms = .True.
   if (body) use_waveforms = .True.
   if (surf) use_waveforms = .True.
   if (dart) use_waveforms = .True.
   end subroutine check_waveforms


   subroutine check_static(static, insar, use_static)
!
!  Args:
!  static: True if static GNSS data used in modelling, False otherwise
!  insar: True if insar data used in modelling, False otherwise
!  use_static: True if some static data used in modelling, False otherwise
!
   implicit none
   logical :: static, insar, use_static
   use_static = .False.
   if (static) use_static = .True.
   if (insar) use_static = .True.
   end subroutine check_static


   subroutine waveform_ffm(strong, cgnss, body, surf, dart, &
       & slip, rake, rupt_time, t_rise, t_fall, many_events)
!
!  Args:
!  strong: True if strong motion data are used, False otherwise
!  cgnss: True if cGNSS data are used, False otherwise
!  body: True if body wave data are used, False otherwise
!  surf: True if surface wave data are used, False otherwise
!  dart: True if DART data are used, False otherwise
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  t_rise: array with model risetime values for all subfaults
!  t_fall: array with model falltime values for all subfaults
!  many_events: True if more than one event to be modelled, False otherwise
!
   implicit none
   real :: slip(:), rake(:), rupt_time(:)
   real :: t_rise(:), t_fall(:)
   logical :: strong, cgnss, dart, body, surf
   logical :: get_coeff, static, insar
   logical :: use_waveforms, many_events
   use_waveforms = .True.
   get_coeff = .True.
   static = .False.
   insar = .False.
   call annealing_set_procedure_param()
   call retrievegf_set_fault_parameters()
   call annealing_set_fault_parameters()
   call saveforward_set_fault_parameters()
   call fourier_coefs()
   call meyer_yamada()
   call get_data(strong, cgnss, body, surf, dart)
   call wavelets_set_data_properties()
   call misfit_eval_set_data_properties()
   call retrievegf_set_data_properties()
   call wavelets_set_data_properties()
   call saveforward_set_data_properties()
   call annealing_set_data_properties()
   call allocate_forward()
   call get_source_fun()
   call get_gf(strong, cgnss, body, surf, dart, many_events)
   call initial_model(slip, rake, rupt_time, t_rise, t_fall)
   t = t_mid
   if (start_annealing .eq. 0) t = t0
   if (many_events) then
      call annealing_set_events()
      call print_summary2(slip, rake, rupt_time, t_rise, t_fall, static, &
        &  insar, get_coeff)
   else
      call print_summary(slip, rake, rupt_time, t_rise, t_fall, static, &
        &  insar, get_coeff)
   endif
   write(*,*)'Start simmulated annealing...'
   do i = 1, n_iter
      if (many_events) then
         call annealing_iter5(slip, rake, rupt_time, t_rise, t_fall, t)
      else
         call annealing_iter3(slip, rake, rupt_time, t_rise, t_fall, t)
      endif
      write(*,*)'iter: ', i
      if (t .lt. t_stop) then
         t = t*0.995
      else
         t = t*cooling_rate
      end if
   end do
   get_coeff = .False.
   if (many_events) then
      call print_summary2(slip, rake, rupt_time, t_rise, t_fall, static, &
        &  insar, get_coeff)
   else
      call print_summary(slip, rake, rupt_time, t_rise, t_fall, static, &
        &  insar, get_coeff)
   endif
   call write_forward(slip, rake, rupt_time, t_rise, t_fall, &
           &  strong, cgnss, body, surf, dart)
   call write_model(slip, rake, rupt_time, t_rise, t_fall, use_waveforms)
   call deallocate_source()
   call deallocate_forward()
   call deallocate_gf()
   end subroutine waveform_ffm
   

   subroutine mixed_ffm(strong, cgnss, body, surf, dart, &
       & static, insar, slip, rake, rupt_time, t_rise, t_fall, many_events)
!
!  Args:
!  strong: True if strong motion data are used, False otherwise
!  cgnss: True if cGNSS data are used, False otherwise
!  body: True if body wave data are used, False otherwise
!  surf: True if surface wave data are used, False otherwise
!  dart: True if DART data are used, False otherwise
!  static: True if static GNSS data used in modelling, False otherwise
!  insar: True if insar data used in modelling, False otherwise
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  rupt_time: array with model rupture time values for all subfaults
!  t_rise: array with model risetime values for all subfaults
!  t_fall: array with model falltime values for all subfaults
!  many_events: True if more than one event to be modelled, False otherwise
!
   implicit none
   real :: slip(:), rake(:), rupt_time(:)
   real :: t_rise(:), t_fall(:)
   logical :: static, strong, cgnss, dart, body, surf
   logical :: get_coeff, insar, ramp_gf_file
   logical :: use_waveforms, many_events
   use_waveforms = .True.
   ramp_gf_file = .False.
   get_coeff = .True.
   call annealing_set_procedure_param()
   call retrievegf_set_fault_parameters()
   call annealing_set_fault_parameters()
   call annealingstatic_set_fault_properties()
   call saveforward_set_fault_parameters()
   call staticdata_set_fault_parameters()
   call insardata_set_fault_parameters()
   call fourier_coefs()
   call meyer_yamada()
   call get_data(strong, cgnss, body, surf, dart)
   call wavelets_set_data_properties()
   call misfit_eval_set_data_properties()
   call retrievegf_set_data_properties()
   call wavelets_set_data_properties()
   call saveforward_set_data_properties()
   call annealing_set_data_properties()
   call allocate_forward()
   call get_source_fun()
   call get_gf(strong, cgnss, body, surf, dart, many_events)
   call initial_model(slip, rake, rupt_time, t_rise, t_fall)
   t = t_mid
   if (start_annealing .eq. 0) t = t0
   if (static) call initial_gnss(slip, rake, many_events)
   if (insar) call get_insar_gf()
   if (insar) call get_insar_data()
   if (insar) call is_ramp(ramp_gf_file)
   if (ramp_gf_file .eqv. .False.) then
      if (insar) call initial_insar(slip, rake)
      if (many_events) then
         call annealing_set_events()
         call print_summary2(slip, rake, rupt_time, t_rise, t_fall, static, &
          &  insar, get_coeff)
      else
         call print_summary(slip, rake, rupt_time, t_rise, t_fall, static, &
          &  insar, get_coeff)
      endif
      write(*,*)'Start simmulated annealing...'
      do i = 1, n_iter
         if (many_events) then
            call annealing_iter6(slip, rake, rupt_time, t_rise, &
            & t_fall, t, static, insar)
         else
            call annealing_iter4(slip, rake, rupt_time, t_rise, t_fall, &
            &  t, static, insar)
         endif
         write(*,*)'iter: ', i
         if (t .lt. t_stop) then
            t = t*0.995
         else
            t = t*cooling_rate
         end if
      end do
      get_coeff = .False.
      if (many_events) then
         call print_summary2(slip, rake, rupt_time, t_rise, t_fall, static, &
          &  insar, get_coeff)
      else
         call print_summary(slip, rake, rupt_time, t_rise, t_fall, static, &
          &  insar, get_coeff)
      endif
      if (insar) call initial_insar(slip, rake)
   else
      call initial_ramp(ramp) 
      call initial_insar(slip, rake, ramp)
      if (many_events) then
         call print_summary2(slip, rake, rupt_time, t_rise, t_fall, static, &
          &  insar, get_coeff, ramp)
      else
         call print_summary(slip, rake, rupt_time, t_rise, t_fall, static, &
          &  insar, get_coeff, ramp)
      endif
      write(*,*)'Start simmulated annealing...'
      do i = 1, n_iter
         if (many_events) then
            call annealing_iter6(slip, rake, rupt_time, t_rise, t_fall, &
            &  t, static, insar, ramp)
         else
            call annealing_iter4(slip, rake, rupt_time, t_rise, t_fall, &
            &  t, static, insar, ramp)
         endif
         write(*,*)'iter: ', i
         if (t .lt. t_stop) then
            t = t*0.995
         else
            t = t*cooling_rate
         end if
      end do
      get_coeff = .False.
      if (many_events) then
         call print_summary2(slip, rake, rupt_time, t_rise, t_fall, static, &
          &  insar, get_coeff, ramp)
      else
         call print_summary(slip, rake, rupt_time, t_rise, t_fall, static, &
          &  insar, get_coeff, ramp)
      endif
      call initial_insar(slip, rake, ramp)
   end if
   call write_forward(slip, rake, rupt_time, t_rise, t_fall, &
        &  strong, cgnss, body, surf, dart)
   if (static) call initial_gnss(slip, rake, many_events)
   call write_model(slip, rake, rupt_time, t_rise, t_fall, use_waveforms)
   call deallocate_source()
   call deallocate_forward()
   call deallocate_gf()
   if (insar) call deallocate_insar_gf()
   end subroutine mixed_ffm


   subroutine static_ffm(slip, rake, static, insar, many_events)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  static: True if static GNSS data used in modelling, False otherwise
!  insar: True if insar data used in modelling, False otherwise
!  many_events: True if more than one event to be modelled, False otherwise
!
   implicit none
   real :: slip(:), rake(:)
   logical :: static, get_coeff, insar, ramp_gf_file
   logical :: use_waveforms, many_events
   use_waveforms = .False.
   get_coeff = .True.
   ramp_gf_file = .False.
   call annealingstatic_set_procedure_param()
   call annealing_set_fault_parameters()
   call annealingstatic_set_fault_properties()
   call initial_model(slip, rake, rupt_time0, t_rise0, t_fall0)
   call staticdata_set_fault_parameters()
   call insardata_set_fault_parameters()
   t = t_mid
   if (start_annealing .eq. 0) t = t0
   if (static) call initial_gnss(slip, rake, many_events)
   if (insar) call get_insar_gf()
   if (insar) call get_insar_data()
   if (insar) call is_ramp(ramp_gf_file)
   if (ramp_gf_file .eqv. .False.) then
      if (insar) call initial_insar(slip, rake)
      call print_static_summary(slip, rake, static, insar, get_coeff)
      write(*,*)'Start simmulated annealing...'
      do i = 1, n_iter
         call annealing_iter(slip, rake, t, static, insar)
         write(*,*)'iter: ', i
         if (t .lt. t_stop) then
            t = t*0.995
         else
            t = t*cooling_rate
         end if
      end do
      get_coeff = .False.
      call print_static_summary(slip, rake, static, insar, get_coeff)
      if (insar) call initial_insar(slip, rake)
   else
      call initial_ramp(ramp) 
      call initial_insar(slip, rake, ramp)
      call print_static_summary(slip, rake, static, insar, get_coeff, ramp)
      write(*,*)'Start simmulated annealing...'
      do i = 1, n_iter
         call annealing_iter(slip, rake, t, static, insar, ramp)
         write(*,*)'iter: ', i
         if (t .lt. t_stop) then
            t = t*0.995
         else
            t = t*cooling_rate
         end if
      end do
      get_coeff = .False.
      call print_static_summary(slip, rake, static, insar, get_coeff, ramp)
      call initial_insar(slip, rake, ramp)
   end if
   if (static) call initial_gnss(slip, rake, many_events)
   call write_model(slip, rake, rupt_time0, t_rise0, t_fall0, use_waveforms)
   if (insar) call deallocate_insar_gf()
   end subroutine static_ffm


end module ffm_methods

