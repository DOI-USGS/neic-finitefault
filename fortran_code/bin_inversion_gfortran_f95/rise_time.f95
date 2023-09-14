module rise_time


   use constants, only : wave_pts2, max_stations, max_rise_time_range, pi, twopi
   use model_parameters, only : get_rise_time
   use wavelet_param, only : get_data_param
   use get_stations_data, only : dt_channel, channels
   implicit none
   complex, allocatable :: source(:, :, :, :)
   integer, parameter, private :: double = kind(1d0)


contains


   subroutine realtr(real1, imag1, n)
!
!  Args:
!  real1: real part of data vector
!  imag1: imaginary part of data vector
!  n: length of data vector
!
   implicit none
   real :: real1(*), imag1(*)
   integer :: n
   integer :: i, i1, i2, half_length, length1, length2
   half_length = 2**(n-1)
   length1 = half_length-1
   length2 = half_length+1
   do i = 1, length1
      i1 = length2+i
      i2 = length2-i
      real1(i1) = real1(i2)
      imag1(i1) = -imag1(i2)
   end do
   imag1(length2) = 0.0
   end subroutine realtr
!

   subroutine fft(real1, imag1, n, sn)
!
!  Args:
!  real1: real part of data vector
!  imag1: imaginary part of data vector
!  n: length of data vector
!  sn: 1 or -1, whether to compute FFT or IFFT
!
   implicit none
   real :: real1(*), imag1(*), sn
   integer :: m(25), n
   real :: norm, imag3, real3, v, imag2, real2
   integer :: i, start, j, j1, k, l, length2, length3, length1, chunks
   length1 = 2**n
   norm = length1
   do i = 1, n
      m(i) = 2**(n-i)
   end do
   do l = 1, n
      chunks = 2**(l-1)
      length2 = length1/chunks
      length3 = length2/2
      k = 0
      do i = 1, chunks
         v = sn*twopi*k/norm
         real2 = cos(v)
         imag2 = sin(v)
         start = length2*(i-1)
         do j = start + 1, start + length3
            j1 = j+length3
            real3 = real1(j1)*real2-imag1(j1)*imag2
            imag3 = real1(j1)*imag2+imag1(j1)*real2
            real1(j1) = real1(j)-real3
            imag1(j1) = imag1(j)-imag3
            real1(j) = real1(j)+real3
            imag1(j) = imag1(j)+imag3
         end do
         do j = 2, n
            if (k.LT.m(j)) exit
            k = k-m(j)
         end do
         k = k+m(j)
      end do
   end do
   k = 0
   do j = 1, length1
      if (k .lt. j) goto 7
      real2 = real1(j)
      imag2 = imag1(j)
      j1 = k+1
      real1(j) = real1(j1)
      imag1(j) = imag1(j1)
      real1(j1) = real2
      imag1(j1) = imag2
 7    do i = 1, n
         if (k .lt. m(i)) exit
         k = k-m(i)
      end do
      k = k+m(i)
   end do
   if (sn .ge. 0.) then
      do i = 1, length1
         real1(i) = real1(i)/norm
         imag1(i) = imag1(i)/norm
      end do
   end if
   end subroutine fft


   subroutine fourier_asym_cosine(omega, t1, t2, source)
!
!  Args:
!  omega: frequency value
!  t1: first parameter of asymmetric cosine function
!  t2: second parameter of asymmetric cosine function
!  source: FFT of asymetric cosine evaluated at frequency omega
!
! analitic fourier transform of asymmetric cosine
! 
   implicit none
   complex source
   real*8 omega, t1, t2
   complex*16 first, second, third, fourth, fifth, z0
      
   z0 = cmplx(1.d0, 0.d0, double)
   if (omega .lt. 1.e-6) then
      first = cmplx(t1+t2,0.d0,double)
   else
      first = cmplx(0.d0, -twopi*omega*(t1 + t2), double)
      first = (exp(first) - z0) * cmplx(0.d0, 0.5d0/pi/omega, double)
   end if
   second = cmplx(0.d0, -pi*(2.d0*omega*t1 + 1.d0), double)
   second = (exp(second) - z0) * cmplx(0.d0, -1.d0/pi/(2.d0*omega + 1.d0/t1), double)
   if (abs(2*t1*omega - 1) .lt. 1.e-6) then
      third = cmplx(-t1, 0.d0, double)
   else
      third = cmplx(0.d0, pi*(1.d0 - 2.d0*omega*t1), double)
      third = (exp(third) - z0) * cmplx(0.d0, 1.d0/pi/(1.d0/t1 - 2.d0*omega), double)
   end if
   if (abs(2*t2*omega - 1) .lt. 1.e-6) then
      fourth = cmplx(t2, 0.d0, double)
   else
      fourth = cmplx(0.d0, pi*(1.d0 - 2.d0*omega*t2), double)
      fourth = (exp(fourth) - z0) * cmplx(0.d0, -1.d0/pi/(1.d0/t2 - 2.d0*omega), double)
   end if
   fifth = cmplx(0.d0, -pi*(2.d0*omega*t2 + 1.d0), double)
   fifth = (exp(fifth) - z0) * cmplx(0.d0, 1.d0/pi/(2.d0*omega + 1.d0/t2), double)

   source = cmplx(0.,-twopi*omega*t1)
   source = exp(source)*(fifth + fourth)
   source = source + third + second
   source = source * cmplx(0.5, 0.) + first
   source = source * cmplx(1./(t1 + t2), 0.d0)  
   end subroutine fourier_asym_cosine


   subroutine get_source_fun()
   implicit none
   real :: dt
   real*8 :: df, t1, t2
   integer :: i, channel, isl, isr, jf
   real :: rise_param1, rise_param2
   integer :: windows, lnpt, jmin, jmax, nlen, max_freq
   allocate(source(wave_pts2, max_stations, max_rise_time_range, max_rise_time_range))
   call get_rise_time(rise_param1, rise_param2, windows)
   call get_data_param(lnpt, jmin, jmax, nlen, max_freq)
!       
! Here, we load into memory, the Fourier transform of rise time function 
!
   jf = 2**(lnpt-1)+1
   do channel = 1, max_stations
      dt = dt_channel(channel)
      if (dt .lt. 1.e-4) cycle
      df = 1.d0/(2**lnpt)/dt
      if (abs(dt - 60.0) .gt. 1.e-4) then
         do isr = 1, windows
            do isl = 1, windows
               t1 = rise_param1+(isl-1)*rise_param2
               t2 = rise_param1+(isr-1)*rise_param2
               t1 = max(dt, t1)
               t2 = max(dt, t2)
               do i = 1, jf
                  call fourier_asym_cosine((i-1)*df, t1, t2, source(i, channel, isl, isr))
               end do
            end do
         end do
      else
         do isr = 1, windows
            do isl = 1, windows
               do i = 1, jf
                  source(i, channel, isl, isr) = cmplx(1.0, 0.0)
               end do
            end do
         end do
      end if
   end do
   end subroutine get_source_fun


   subroutine deallocate_source()
   deallocate(source)
   end subroutine deallocate_source


end module rise_time
