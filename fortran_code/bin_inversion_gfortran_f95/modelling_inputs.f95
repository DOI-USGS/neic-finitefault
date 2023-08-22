module modelling_inputs


   implicit none
   integer :: seed0, io_data, n_iter
   real :: t0, cooling_rate, t_stop, emin0, smooth_moment, smooth_slip, smooth_time
   real :: t_mid, t_latest, moment_input
   integer :: io_re, io_func, io, events
   real :: moment_event(10)


contains


   subroutine read_annealing_param()
   implicit none
   write(*,*)'Store annealing and modelling parameters...'
   open(10, file='annealing.txt', status='old')
   read(10,*) n_iter, seed0, io_data, moment_input
   read(10,*) t0, cooling_rate, t_stop, emin0, smooth_moment, smooth_slip, smooth_time
   read(10,*) io_re, t_mid, io_func, t_latest
   read(10,*) io
   close(10)
   end subroutine read_annealing_param
   
   
   subroutine moment_events()
   implicit none
   integer :: io, i
   open(10, file='moment_events.txt', status='old')
   events = 0
   do
      read(10,*,iostat=io)
      if (io/=0) exit
      events = events + 1
   enddo
   rewind (10)
   do i=1, events
      read(10,*)moment_event(i)
   enddo
   close(10)
   end subroutine moment_events
   

   subroutine get_annealing_param(n_iter0, seed00, t00, cooling_rate0, t_stop0, io_re0, t_mid0)
!
!  Args:
!  n_iter0: iterations of annealing method
!  seed00: initial seed for random number generators
!  t00: initial temperature for annealing method
!  cooling_rate0: cooling rate of annealing method
!  t_stop0: temperature at which cooling of annealing method becomes slower
!  io_re0: whether to perform simmulated annealing
!  t_mid0: another temperature value
!  
   implicit none
   integer :: n_iter0, seed00, io_re0
   real :: t00, cooling_rate0, t_stop0, t_mid0
   n_iter0 = n_iter
   seed00 = seed0
   t00 = t0
   cooling_rate0 = cooling_rate
   t_stop0 = t_stop
   io_re0 = io_re
   t_mid0 = t_mid
   end subroutine get_annealing_param
   

   subroutine get_weights_moment_end(moment_input0, smooth_moment0, smooth_slip0, smooth_time0, &
                                &    t_latest0, emin00)
!
!  Args:
!  moment_input0: input seismic moment
!  smooth_moment0: weight of moment constrain
!  smooth_slip0: weight of slip constrain
!  smooth_time0: weight of time constrain
!  t_latest0: latest time of rupture time allowed
!  emin00: related to moment constrain
!  
   implicit none
   real :: moment_input0, smooth_moment0, smooth_slip0, smooth_time0, t_latest0, emin00
   moment_input0 = moment_input
   smooth_moment0 = smooth_moment
   smooth_slip0 = smooth_slip
   smooth_time0 = smooth_time
   t_latest0 = t_latest
   emin00 = emin0
   end subroutine get_weights_moment_end
   

   subroutine get_many_events(events0, moment_event0)
!
!  Args:
!  events0: amount of events to model
!  moment_event0: seismic moment for each event
!  
   implicit none
   integer :: events0
   real :: moment_event0(10)
   events0 = events
   moment_event0(:) = moment_event(:)
   end subroutine get_many_events
   

end module modelling_inputs 
