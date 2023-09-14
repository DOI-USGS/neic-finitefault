program run_modelling


   use constants, only : max_seg, max_subfaults, max_subf
   use model_parameters, only : get_faults_data, get_model_space, get_special_boundaries, &
         & subfault_positions, modelparameters_set_procedure_param, &
         & modelparameters_set_events, events_segments, deallocate_ps
   use modelling_inputs, only : read_annealing_param, moment_events, get_annealing_param
   use save_forward, only : write_forward
   use random_gen, only : start_seed
   use annealing, only : n_threads
   use ffm_methods, only: check_waveforms, check_static, waveform_ffm, &
                     &  mixed_ffm, static_ffm, ffmmethod_set_procedure_param
   use regularization, only : regularization_set_fault_parameters 
   implicit none
   integer :: i
   real :: slip(max_subfaults), rake(max_subfaults), rupt_time(max_subfaults)
   real :: t_rise(max_subfaults), t_fall(max_subfaults)
   real :: t
   real*8 :: ramp(18)
   logical :: static, strong, cgps, dart, body, surf, auto
   logical :: insar, ramp_gf_file
   logical :: use_waveforms, use_static, many_events
   integer :: int0, seed0, int1
   real :: real0, real1, real2, real3
   character(len=10) :: input

   write(*,'(/A/)')"CHEN-JI'S WAVELET KINEMATIC MODELLING METHOD"
   static = .False.
   insar = .False.
   strong = .False.
   cgps = .False.
   body = .False.
   surf = .False.
   dart = .False.
   auto = .False.
   many_events = .False.
   do i = 1, iargc()
      call getarg(i, input)
      input = trim(input)
      if (input .eq.'gps') static = .True.
      if (input .eq.'insar') insar = .True.
      if (input .eq.'strong') strong = .True.
      if (input .eq.'cgps') cgps = .True.
      if (input .eq.'body') body = .True.
      if (input .eq.'surf') surf = .True.
      if (input .eq.'dart') dart = .True.
      if (input .eq.'auto') auto = .True.
      if (input .eq.'many') many_events = .True.
   end do
   call check_waveforms(strong, cgps, body, surf, dart, use_waveforms)
   call check_static(static, insar, use_static)
   call n_threads(auto)
   call read_annealing_param()
   if (many_events) call moment_events()
   if (many_events) call modelparameters_set_events()
   call get_annealing_param(int0, seed0, real0, real1, real2, int1, real3)
   call start_seed(seed0)
   call modelparameters_set_procedure_param()
   call get_faults_data()
   call get_model_space()
   call get_special_boundaries()
   call subfault_positions()
   if (many_events) call events_segments()
   call regularization_set_fault_parameters()
   call ffmmethod_set_procedure_param()
   if ((use_waveforms) .and. (use_static .eqv. .False.)) then
      call waveform_ffm(strong, cgps, body, surf, dart, &
       & slip, rake, rupt_time, t_rise, t_fall, many_events)
   elseif ((use_static) .and. (use_waveforms .eqv. .False.)) then
      call static_ffm(slip, rake, static, insar, many_events)
   elseif ((use_static) .and. (use_waveforms)) then
      call mixed_ffm(strong, cgps, body, surf, dart, &
       & static, insar, slip, rake, rupt_time, t_rise, t_fall, many_events)
   endif
   call deallocate_ps()
   write(*,'(/A/)')"END CHEN-JI'S WAVELET KINEMATIC MODELLING METHOD"


end program run_modelling

