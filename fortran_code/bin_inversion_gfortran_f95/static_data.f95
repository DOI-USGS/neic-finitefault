!     This program load the static subfault GFs and generate the matrix
!     for baseline correction.         Ji, Chen, 2003
!     Reviewed. The previous version have an error in loading GFs.
!                      
!                                   Ji, Chen, 2004
!  Fixed
!
!
!
module static_data


   use constants, only : max_subf, max_seg, max_subfaults, dpi
   implicit none
   integer, parameter, private :: n_stations = 500
   integer :: n_chan
   real*8 :: synm_whole(n_stations, 3), weight_sum
   real :: lat(n_stations), lon(n_stations), max_gnss
   real :: green(n_stations, 6, max_subfaults), syn_disp(n_stations, 3)
   real :: obse(n_stations, 3), weight(n_stations, 3)
   character(len=6) :: sta_name(n_stations)
   integer :: subfaults, segments, nxs_sub(max_seg), nys_sub(max_seg)
   logical :: subfault_in_event(max_subfaults, 10)
   integer :: event_sta(n_stations)
  

contains


   subroutine staticdata_set_fault_parameters()
   use model_parameters, only : get_segments, get_events_segments
   implicit none
   real :: dip(max_seg), strike(max_seg), delay_seg(max_seg)
   integer :: cum_subfaults(max_seg)
   logical :: segment_in_event(max_seg, 10)
   call get_segments(nxs_sub, nys_sub, dip, strike, delay_seg, segments, subfaults, cum_subfaults)
   call get_events_segments(segment_in_event, subfault_in_event)
   end subroutine staticdata_set_fault_parameters
   

   subroutine initial_gnss(slip, rake, many_events)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  many_events: True if modelling more than a single earthquake, False otherwise
!
   implicit none
   integer k, j, segment, channel, i
   integer iys, ixs, int1, subfault, event
   real slip(:), rake(:)
   real :: cosal, sinal, angle
   real*8 :: disp
   logical is_file, is_file2, many_events
   character(len=30) :: event_file, string1, string2
!
   inquire( file = 'static_data.txt', exist = is_file )
   if (is_file) then 
      weight_sum = 0.d0
      open(9, file='static_data.txt', status='old')
      read(9,*) n_chan
      read(9,*)
      do i = 1, n_chan
         read(9,*) int1, sta_name(i), lat(i), lon(i), (obse(i, k), k = 1, 3), (weight(i, k), k = 1, 3)
         do k = 1, 3
            weight_sum = weight_sum + weight(i, k)
         end do
      end do
      max_gnss = maxval(abs(obse(:, :)))
      close(9)

      event_file = 'static_events.txt'
      inquire(file = event_file, exist = is_file)
      if (is_file) then
         open(12, file=event_file, status='old')
         do i=1,n_chan
            read(12,*)string1, string2, event_sta(i)
         enddo
         close(12)    
      endif

      open(33, file='Green_static_subfault.txt', status='old')
      read(33,*) int1
      do channel = 1, n_chan
         if (many_events) event = event_sta(channel)
         read(33,*)
         subfault = 0
         do segment = 1, segments
            read(33,*)
            do j = 1, nys_sub(segment)*nxs_sub(segment)
               subfault = subfault + 1
               read(33,*)(green(channel, k, subfault), k = 1, 6)
               if (many_events .and. (subfault_in_event(subfault, event) .eqv. .False.)) then
                  green(channel, :, subfault) = 0.0
               endif
            enddo
         end do
      end do
      close(33)
      
      do k = 1, 3
         j = 2 * k - 1
         do channel = 1, n_chan
            disp = 0.d0
            do subfault = 1, subfaults
               angle = rake(subfault)*dpi 
               sinal = sin(angle)
               cosal = cos(angle)
               disp = disp + slip(subfault) &
     &  *(sinal*green(channel, j, subfault)+cosal*green(channel, j+1, subfault))
            end do
            syn_disp(channel, k) = disp
         end do
      end do 
      open(10,file='static_synthetics.txt')
      write(10,*) n_chan
      do channel = 1, n_chan
         write(10,*) channel, sta_name(channel), lat(channel), lon(channel),(syn_disp(channel, k), k=1, 3)
      end do
      close(10)
   end if
   end subroutine initial_gnss
   
   
!
! routine for loading static synthetic seismograms, given a rupture model
!
   subroutine static_synthetic(slip, rake, err)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  err: static GNSS misfit
!
   implicit none
   real slip(:), rake(:)!, err
   integer k, j, segment, channel, subfault
   real err, dif, angle, sinal, cosal
   real*8 :: disp, err2
   
   err2 = 0.d0
   synm_whole(:, :) = 0.d0
   do k = 1, 3
      j = 2 * k - 1
      do channel = 1, n_chan
         disp = 0.d0
         do subfault = 1, subfaults
            angle = rake(subfault)*dpi 
            sinal = sin(angle)
            cosal = cos(angle)
            disp = disp + slip(subfault) &
       &  *(sinal*green(channel, j, subfault)+cosal*green(channel, j+1, subfault))
         end do
         synm_whole(channel, k) = disp
         dif = synm_whole(channel, k) - obse(channel, k)
         err2 = err2 + weight(channel, k) * dif * dif / max_gnss / max_gnss!100.0
      end do
   end do
   err2 = sqrt((err2/weight_sum))
   err = real(err2)
   end subroutine static_synthetic
            

!  
! subroutine for removing the static response of current subfault for all static stations
!  
   subroutine static_remove_subfault(slip, rake, subfault)
!
!  Args:
!  slip: slip value for this subfault
!  rake: rake value for this subfault
!  subfault: value of current subfault
!
   implicit none
   real, intent(in) :: slip, rake
   integer, intent(in) :: subfault
   integer k, j, channel
   real disp, angle, sinal, cosal
!
   angle = rake*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   do k = 1, 3
      j = 2*k-1
      do channel = 1, n_chan
         disp = slip &
     &   *(sinal*green(channel, j, subfault)+cosal*green(channel, j+1, subfault))
         synm_whole(channel, k) = synm_whole(channel, k)-disp
      end do
   end do

   end subroutine static_remove_subfault


!
! subroutine for testing the new response of the current subfault, for all stations
! we also give the misfit error of static data
!
   pure subroutine static_modify_subfault(slip, rake, subfault, err)
!
!  Args:
!  slip: slip value for this subfault
!  rake: rake value for this subfault
!  subfault: value of current subfault
!  err: static GNSS misfit
!
   implicit none
   real, intent(in) :: slip, rake
   real, intent(out) :: err
   integer, intent(in) :: subfault
   integer k, j, channel
   real disp, angle, dif, sinal, cosal
   real*8 :: err2
!
   angle = rake*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   err2 = 0.d0
   do k = 1, 3
      j = 2*k-1
      do channel = 1, n_chan
         disp = slip &
       & *(sinal*green(channel, j, subfault)+cosal*green(channel, j+1, subfault))
         dif = (synm_whole(channel, k) + disp) - obse(channel, k)
         err2 = err2 + weight(channel, k) * dif * dif / max_gnss / max_gnss!100.0
      end do
   end do
   err2 = sqrt(err2/weight_sum)
   err = real(err2)

   end subroutine static_modify_subfault


!
! subroutine for asliping response of current subfault, for all stations.
! we also give the misfit error of static data
!
   subroutine static_add_subfault(slip, rake, subfault)
!
!  Args:
!  slip: slip value for this subfault
!  rake: rake value for this subfault
!  subfault: value of current subfault
!
   implicit none
   real, intent(in) :: slip, rake
   integer, intent(in) :: subfault
   integer k, j, channel
   real disp, angle, dif, sinal, cosal
   
   angle = rake*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   do k = 1, 3
      j = 2*k-1
      do channel = 1, n_chan
         disp = slip &
       & *(sinal*green(channel, j, subfault)+cosal*green(channel, j+1, subfault))
         synm_whole(channel, k) = synm_whole(channel, k)+disp
      end do
   end do

   end subroutine static_add_subfault


end module static_data
