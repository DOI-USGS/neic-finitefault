module rise_time


   use constants, only : wave_pts2, max_stations, max_rise_time_range, pi, twopi
   use model_parameters, only : get_rise_time
   use wavelet_param, only : get_data_param
   use get_stations_data, only : dt_channel, channels
   implicit none
   complex, allocatable :: source(:, :, :, :)
   integer, parameter, private :: double = kind(1d0)


contains


   subroutine realtr(xr, xi, n)
!
!  Args:
!  xr: real part of data vector
!  xi: imaginary part of data vector
!  n: length of data vector
!
   implicit none
   real :: xr(*), xi(*)
   integer :: n
   integer :: i, i1, i2, lb, lh
   lh = 2**(n-1)
   lb = lh-1
   lh = lh+1
   do i = 1, lb
      i1 = lh+i
      i2 = lh-i
      xr(i1) = xr(i2)
      xi(i1) = -xi(i2)
   end do
   xi(lh) = 0.0
   end subroutine realtr
!

   subroutine fft(xr, xi, n, sn)
!
!  Args:
!  xr: real part of data vector
!  xi: imaginary part of data vector
!  n: length of data vector
!  sn: 1 or -1, whether to compute FFT or IFFT
!
   implicit none
   real :: xr(*), xi(*), sn
   integer :: m(25), n
   real :: flx, holdi, holdr, qi, qr, v, wki, wkr
   integer :: i, ib, ist, j, j1, jh, k, l, lb, lbh, lx, nb 
   lx = 2**n
   flx = lx
   do i = 1, n
      m(i) = 2**(n-i)
   end do
   do l = 1, n
      nb = 2**(l-1)
      lb = lx/nb
      lbh = lb/2
      k = 0
      do ib = 1, nb
         v = sn*twopi*k/flx
         wkr = cos(v)
         wki = sin(v)
         ist = lb*(ib-1)
         do i = 1, lbh
            j = ist+i
            jh = j+lbh
            qr = xr(jh)*wkr-xi(jh)*wki
            qi = xr(jh)*wki+xi(jh)*wkr
            xr(jh) = xr(j)-qr
            xi(jh) = xi(j)-qi
            xr(j) = xr(j)+qr
            xi(j) = xi(j)+qi
         end do
         do i = 2, n
            if (k.LT.m(i)) exit
            k = k-m(i)
         end do
         k = k+m(i)
      end do
   end do
   k = 0
   do j = 1, lx
      if (k .lt. j) goto 7
      holdr = xr(j)
      holdi = xi(j)
      j1 = k+1
      xr(j) = xr(j1)
      xi(j) = xi(j1)
      xr(j1) = holdr
      xi(j1) = holdi
 7    do i = 1, n
         if (k .lt. m(i)) exit
         k = k-m(i)
      end do
      k = k+m(i)
   end do
   if (sn .ge. 0.) then
      do i = 1, lx
         xr(i) = xr(i)/flx
         xi(i) = xi(i)/flx
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
   integer :: i, ir, isl, isr, jf
   real :: ta0, dta
   integer :: msou, lnpt, jmin, jmax, nlen, max_freq
   allocate(source(wave_pts2, max_stations, max_rise_time_range, max_rise_time_range))
   call get_rise_time(ta0, dta, msou)
   call get_data_param(lnpt, jmin, jmax, nlen, max_freq)
!       
! Here, we load into memory, the Fourier transform of rise time function 
!
   jf = 2**(lnpt-1)+1
   do ir = 1, max_stations
      dt = dt_channel(ir)
      if (dt .lt. 1.e-4) cycle
      df = 1.d0/(2**lnpt)/dt
      if (abs(dt - 60.0) .gt. 1.e-4) then
         do isr = 1, msou
            do isl = 1, msou
               t1 = ta0+(isl-1)*dta
               t2 = ta0+(isr-1)*dta
               t1 = max(dt, t1)
               t2 = max(dt, t2)
               do i = 1, jf
                  call fourier_asym_cosine((i-1)*df, t1, t2, source(i, ir, isl, isr))
               end do
            end do
         end do
      else
         do isr = 1, msou
            do isl = 1, msou
               do i = 1, jf
                  source(i, ir, isl, isr) = cmplx(1.0, 0.0)
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
