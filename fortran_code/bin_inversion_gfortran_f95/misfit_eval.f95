module misfit_eval


   use constants, only : wave_pts2, max_stations
   use wavelet_param, only : get_data_param
   use get_stations_data, only : get_options, get_wavelet_obs
   implicit none
   real :: weight(max_stations), wave_obs(wave_pts2, max_stations), wmax(max_stations)
   integer :: misfit_type(12, max_stations), lnpt, jmin, jmax, nlen
   integer :: t_min(max_stations), t_max(max_stations)
   real :: wavelet_weight(12, max_stations)


contains


   subroutine misfit_eval_set_data_properties()
   implicit none
   integer :: max_freq, t_max_val(max_stations)
   call get_data_param(lnpt, jmin, jmax, nlen, max_freq)
   call get_options(weight, misfit_type, t_min, t_max, wavelet_weight)
   call get_wavelet_obs(wave_obs, wmax, t_max_val)
   end subroutine misfit_eval_set_data_properties


   pure subroutine misfit_channel(channel, wave_syn, error2)
!
!  Args:
!  channel: number of current channel
!  wave_syn: Wavelet coefficients of synthetic waveform for current channel
!  error2: misfit between synthetic and observed wavelet coeficients for this channel
!  
!  Misfit between observed and synthetic waveforms, in wavelet domain.
!
   real, intent(inout) :: wave_syn(wave_pts2)
   integer, intent(in) :: channel
   real*8, intent(out) :: error2
   real*8 :: aa, ab, bb, er, er1, er2, erp!, error2
   real :: ramp
   integer :: i, j, k, n1, n2, n_be, n_begin, n_delt, n_ed
!     
   error2 = 0.d0
   ramp = wmax(channel)
   n1 = 2**(jmin-1)
   n2 = 2**(jmax)-1
   do i = n1, n2
      wave_syn(i) = wave_syn(i)/ramp
   end do
!   ramp = 1
   do j = jmin, jmax
      er = 0.d0
      if (misfit_type(j, channel) .eq. 0) cycle
      if (wavelet_weight(j, channel) .lt. 1.0e-5) cycle
      n_begin = 2**(j-1)
      n_delt = nlen/n_begin
      n_be = n_begin+int(t_min(channel)/n_delt)
      n_ed = n_begin+int(t_max(channel)/n_delt+0.5)-1
      if (n_ed .lt. n_be) n_ed = n_be
!  j = 1 L1 Norm
      if (misfit_type(j, channel) .eq. 1) then
         do k = n_be, n_ed
            er = er+abs(wave_syn(k)-wave_obs(k, channel))
         end do
         er = wavelet_weight(j, channel)*(er/(n_ed-n_be+1))
         error2 = error2 + er
      end if
!       j = 2 L2 norm
      if (misfit_type(j, channel) .eq. 2) then
         do k = n_be, n_ed
            erp = wave_syn(k)-wave_obs(k, channel)
            er = er+erp*erp
         end do
         er = wavelet_weight(j, channel)*sqrt(er/(n_ed-n_be+1))
         error2 = error2 + er
      end if
!       j = 3 L1+L2 Norm
      if (misfit_type(j, channel) .eq. 3) then
         er1 = 0.d0
         er2 = 0.d0
         do k = n_be, n_ed
            erp = abs(wave_syn(k)-wave_obs(k, channel))
            er1 = er1+erp 
            er2 = er2+erp*erp
         end do
         er = wavelet_weight(j, channel) &
      &  *(0.5d0*sqrt(er2 / (n_ed-n_be+1))+0.5d0*er1 / (n_ed-n_be+1))
         error2 = error2 + er
      end if
!     j = 4 correlation
      if (misfit_type(j, channel) .eq. 4) then
         ab = 0.d0
         aa = 0.d0
         bb = 0.d0
         do k = n_be, n_ed
            ab = wave_syn(k)*wave_obs(k, channel)+ab
            aa = wave_syn(k)*wave_syn(k)+aa
            bb = wave_obs(k, channel)*wave_obs(k, channel)+bb
         end do
         er = wavelet_weight(j, channel)*(1.d0-2.d0*ab/(aa+bb))
         error2 = error2 + er
      end if
   end do
   error2 = error2/(jmax-jmin+1)
   error2 = error2*weight(channel)
   end subroutine misfit_channel


end module misfit_eval
