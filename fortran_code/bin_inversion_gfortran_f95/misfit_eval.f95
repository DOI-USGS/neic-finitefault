module misfit_eval


   use constants, only : wave_pts2, max_stations
   use wavelet_param, only : get_data_param
   use get_stations_data, only : get_options, get_wavelet_obs
   implicit none
   real :: weight(max_stations), wave_obs(wave_pts2, max_stations), max_coeff(max_stations)
   integer :: misfit_type(12, max_stations), lnpt, jmin, jmax, nlen
   integer :: t_min(max_stations), t_max(max_stations)
   real :: wavelet_weight(12, max_stations)


contains


   subroutine misfit_eval_set_data_properties()
   implicit none
   integer :: max_freq, index_max(max_stations)
   call get_data_param(lnpt, jmin, jmax, nlen, max_freq)
   call get_options(weight, misfit_type, t_min, t_max, wavelet_weight)
   call get_wavelet_obs(wave_obs, max_coeff, index_max)
   end subroutine misfit_eval_set_data_properties


   pure subroutine misfit_channel(channel, wave_syn, misfit)
!
!  Args:
!  channel: number of current channel
!  wave_syn: Wavelet coefficients of synthetic waveform for current channel
!  misfit: misfit between synthetic and observed wavelet coeficients for this channel
!  
!  Misfit between observed and synthetic waveforms, in wavelet domain.
!
   real, intent(inout) :: wave_syn(wave_pts2)
   integer, intent(in) :: channel
   real*8, intent(out) :: misfit
   real*8 :: misfit0, value1, value2, value3, diff!, misfit 
   real :: max_coeff0
   integer :: i, j, k, n1, n2, start_scale, scale_size
   integer :: misfit_type0, scale_weight
!     
   misfit = 0.d0
   max_coeff0 = max_coeff(channel)
   n1 = 2**(jmin-1)
   n2 = 2**(jmax)-1
   do i = n1, n2
      wave_syn(i) = wave_syn(i)/max_coeff0
   end do
!   max_coeff = 1
   do j = jmin, jmax
      misfit_type0 = misfit_type(j, channel)
      scale_weight = wavelet_weight(j, channel)
      misfit0 = 0.d0
      value1 = 0.d0
      value2 = 0.d0
      value3 = 0.d0
      if (misfit_type0 .eq. 0) cycle
      if (scale_weight .lt. 1.0e-5) cycle
      start_scale = 2**(j-1)
      scale_size = nlen/start_scale
      n1 = start_scale+int(t_min(channel)/scale_size)
      n2 = start_scale+int(t_max(channel)/scale_size+0.5)-1
      if (n2 .lt. n1) n2 = n1
!  j = 1 L1 Norm
      select case (misfit_type0)
         case (1)      
            do k = n1, n2
               value1 = value1+abs(wave_syn(k)-wave_obs(k, channel))
            end do
            misfit0 = scale_weight*(value1/(n2-n1+1))
         case (2)
            do k = n1, n2
               diff = wave_syn(k)-wave_obs(k, channel)
               value2 = value2+diff*diff
            end do
            misfit0 = scale_weight*sqrt(misfit0/(n2-n1+1))
         case (3)
            do k = n1, n2
               diff = abs(wave_syn(k)-wave_obs(k, channel))
               value1 = value1+diff 
               value2 = value2+diff*diff
            end do
            misfit0 = 0.5d0*scale_weight &
         &  *(sqrt(value2 / (n2-n1+1))+value1 / (n2-n1+1))
         case (4)
            do k = n1, n2
               value1 = wave_syn(k)*wave_obs(k, channel)+value1
               value2 = wave_syn(k)*wave_syn(k)+value2
               value3 = wave_obs(k, channel)*wave_obs(k, channel)+value3
            end do
            misfit0 = scale_weight*(1.d0-2.d0*value3/(value1+value2))
      end select
      misfit = misfit + misfit0
   end do
   misfit = misfit/(jmax-jmin+1)
   misfit = misfit*weight(channel)
   end subroutine misfit_channel


end module misfit_eval
