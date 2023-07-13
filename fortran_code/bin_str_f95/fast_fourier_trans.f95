!c
module fast_fourier_trans


   use constants, only : twopi
   implicit none
   real, private :: cos_fft(4200), sin_fft(4200)!, 15), sin_fft(4200, 15)
   integer, private :: kkk(4200,15)


contains


   subroutine fourier_coefs()
   implicit none
!
! in this subroutine, we load into memory certain coeficients and values which are frequently used 
! in the method of the cfft. 
!
   integer :: i, n, ib, k, power, power3, power2
   real*8 :: omega

   cos_fft(:) = 0.
   sin_fft(:) = 0.
   kkk(:, :) = 0
        
   do n =2, 12
      power2 = 2 ** n
      k=0
      power = 2 ** (n - 1)
      do ib = 1, power
         omega = twopi*dble(k)/dble(power2)
         cos_fft(ib) = cos(omega)
         sin_fft(ib) = sin(omega)
         do i = 2, n
            power3 = 2 ** (n - i)
            if(k.lt.power3) exit
            k = k - power3
         enddo
         k = k + power3
      enddo
      k = 0
      do ib = 1, power2
         kkk(ib, n) = k
         do i = 1, n
            power3 = 2 ** (n - i)
            if(k.lt.power3) exit
            k = k - power3
         enddo
         k = k + power3
      enddo
   enddo
   end subroutine fourier_coefs


   subroutine fft(xr, xi, n, sn)
!
!  Args:
!  xr: real part of data vector
!  xi: imaginary part of data vector
!  n: length of data vector
!  sn: 1 or -1, whether to compute FFT or IFFT
!
   implicit none
   integer :: n, sn
   real :: xr(*), xi(*)
   integer :: lx, l, lb, lbh, ib, ist, j, j1, k, nb
   real :: flx, qr, qi, holdr, holdi
   real :: real_exp, imag_exp
   LX = 2 ** N
   FLX = real(lx)
   DO L = 1, N
      NB = 2 ** (L - 1)
      LB = LX / NB
      LBH = LB / 2
      DO IB = 1, NB
         real_exp = cos_fft(ib)
         imag_exp = real(sn) * sin_fft(ib)
         IST = LB * (IB - 1)
         do J = 1 + IST, LBH + IST
            j1 = J + LBH
            QR = XR(j1) * real_exp - XI(j1) * imag_exp
            QI = XR(j1) * imag_exp + XI(j1) * real_exp
            XR(j1) = XR(J) - QR
            XI(j1) = XI(J) - QI
            XR(J) = XR(J) + QR
            XI(J) = XI(J) + QI
         enddo
      enddo
   enddo
   do J = 1, LX
      k = kkk(j, n)
      if(k.lt.j) cycle
      HOLDR = XR(J)
      HOLDI = XI(J)
      J1 = K + 1
      XR(J) = XR(J1)
      XI(J) = XI(J1)
      XR(J1) = HOLDR
      XI(J1) = HOLDI
   enddo
   if(sn.gt.0) then
      xr(:lx) = xr(:lx) / flx
      xi(:lx) = xi(:lx) / flx
   endif
   end subroutine fft


end module fast_fourier_trans
