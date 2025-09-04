program run_modelling


   use constants, only : max_seg, nt1, nnxy
   use model_parameters, only : get_faults_data, get_model_space, get_special_boundaries, subfault_positions, &
                            &   write_model
   use modelling_inputs, only : get_annealing_param, n_iter, io_re, cooling_rate, t_stop, t_mid, t0, io_data, idum
   use get_stations_data, only : get_data
   use retrieve_gf, only : get_gf, deallocate_gf
   use save_forward, only : write_forward
   use rise_time, only : get_source_fun, deallocate_source
   use static_data, only : initial_gnss
   use random_gen, only : start_seed
   use annealing, only : initial_model, initial_regularization, annealing_iter, annealing_iter2, &
                     &   annealing_iter3, annealing_iter4, n_threads
   implicit none
   integer :: i
   real :: slip(nnxy, max_seg), rake(nnxy, max_seg), rupt_time(nnxy, max_seg)
   real :: tl(nnxy, max_seg), tr(nnxy, max_seg)
   real :: er, t
   logical :: static, strong, cgnss, dart, body, surf, auto
   character(len=10) :: input

   static = .False.
   strong = .False.
   cgnss = .False.
   body = .False.
   surf = .False.
   dart = .False.
   auto = .False.
   do i = 1, iargc()
      call getarg(i, input)
      input = trim(input)
      if (input .eq.'gnss') static = .True.
      if (input .eq.'strong') strong = .True.
      if (input .eq.'cgnss') cgnss = .True.
      if (input .eq.'body') body = .True.
      if (input .eq.'surf') surf = .True.
      if (input .eq.'dart') dart = .True.
      if (input .eq.'auto') auto = .True.
   end do
   call n_threads(auto)
   call get_annealing_param()
   call start_seed(idum)
   call get_faults_data()
   call get_model_space()
   call get_special_boundaries()
   call subfault_positions()
   write(*,*) io_data
   call get_gf(strong, cgnss, body, surf, dart)
   call get_data(strong, cgnss, body, surf, dart)
   call get_source_fun()
   call initial_model(slip, rake, rupt_time, tl, tr)
   if (static) call initial_gnss(slip, rake)
   call initial_regularization(slip, rake, rupt_time, tl, tr, static)
   t = t_mid
   if (io_re .eq. 0) t = t0
   if (static) then
      do i = 1, n_iter
         call annealing_iter4(slip, rake, rupt_time, tl, tr, er, t)
         write(*,*)'iter: ', i
         if (t .lt. t_stop) then
            t = t*0.995
         else
            t = t*cooling_rate
         end if
      end do
   else
      do i = 1, n_iter
         call annealing_iter3(slip, rake, rupt_time, tl, tr, er, t)
         write(*,*)'iter: ', i
         if (t .lt. t_stop) then
            t = t*0.995
         else
            t = t*cooling_rate
         end if
      end do
   end if
   call write_forward(slip, rake, rupt_time, tl, tr, strong, cgnss, body, surf, dart)
   call initial_gnss(slip, rake)
   call write_model(slip, rake, rupt_time, tl, tr)
   call deallocate_source()
   call deallocate_gf()


end program run_modelling

