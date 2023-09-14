program run_forward


   use model_parameters, only : get_faults_data, slip0, rake0, rupt_time0, &
       &  t_rise0, t_fall0, write_model, deallocate_ps
   use modelling_inputs, only : read_annealing_param
   use get_stations_data, only : get_data
   use retrieve_gf, only : get_gf, deallocate_gf, retrievegf_set_data_properties, &
                    &   retrievegf_set_fault_parameters
   use wavelets, only : wavelets_set_data_properties, fourier_coefs, meyer_yamada
   use save_forward, only : write_forward, saveforward_set_fault_parameters, &
                    &   saveforward_set_data_properties
   use static_data, only : initial_gps, staticdata_set_fault_parameters
   use insar_data, only : initial_insar, get_insar_gf, deallocate_insar_gf, &
                    &   get_insar_data, is_ramp, initial_ramp, &
                    &   insardata_set_fault_parameters
   implicit none
   integer i
   character(len=10) :: input
   logical :: static, strong, cgps, body, surf, dart, insar
   logical :: use_waveforms, many_events

   static = .False.
   insar = .False.
   strong = .False.
   cgps = .False.
   body = .False.
   surf = .False.
   dart = .False.
   use_waveforms = .True.
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
      if (input .eq.'many') many_events = .True.
   end do
   call read_annealing_param()
   call get_faults_data()
   call retrievegf_set_fault_parameters()
   call saveforward_set_fault_parameters()
   call fourier_coefs()
   call meyer_yamada()
   call get_data(strong, cgps, body, surf, dart)
   call wavelets_set_data_properties()
   call retrievegf_set_data_properties()
   call wavelets_set_data_properties()
   call saveforward_set_data_properties()
   call get_gf(strong, cgps, body, surf, dart, many_events)
   call staticdata_set_fault_parameters()
   call insardata_set_fault_parameters()
   call write_forward(slip0, rake0, rupt_time0, t_rise0, t_fall0, &
       &  strong, cgps, body, surf, dart)
   if (static) call initial_gps(slip0, rake0, many_events)
   if (insar) call get_insar_gf()
   if (insar) call get_insar_data()
   if (insar) call initial_insar(slip0, rake0)
   call write_model(slip0, rake0, rupt_time0, t_rise0, t_fall0, use_waveforms)
   call deallocate_gf()
   call deallocate_ps()
   if (insar) call deallocate_insar_gf()


end program run_forward
