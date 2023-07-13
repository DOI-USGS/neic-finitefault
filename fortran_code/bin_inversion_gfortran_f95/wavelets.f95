module wavelets


   use constants, only : pi, twopi, wave_pts2
   use wavelet_param, only : get_data_param
   implicit none
   complex :: rwt1(wave_pts2, 12), rwt2(wave_pts2, 12), c1, c2
   integer :: kkk(4200, 15)
   real :: cos_fft(4200), sin_fft(4200)
   real :: wavelet_r1(wave_pts2, 12), wavelet_r2(wave_pts2, 12)
   real :: wavelet_i1(wave_pts2, 12), wavelet_i2(wave_pts2, 12)
   real :: real_c1, imag_c1, real_c2, imag_c2
   integer :: jmin, jmax, lnpt, nlen


contains
   

   subroutine wavelets_set_data_properties()
   implicit none
   integer :: max_freq
   call get_data_param(lnpt, jmin, jmax, nlen, max_freq)
   end subroutine wavelets_set_data_properties


   subroutine wavelet_obs(cr, cz, u)
!
!  Args:
!  cr: real part of data vector
!  cz: imaginary part of data vector
!  u: wavelet coefficients from data vector
!
   implicit none
!       ******* Explaination ********
!       In order to normalized the amplitude, we move the coeffficient
!       sqrt(T/kmax)
!                                       Jichen 1997, 1, 20
!
! we split between this routine and wavelet_syn routine to hopefully reduce
! computation time
!
   real :: cr(wave_pts2), cz(wave_pts2), u(wave_pts2)
   integer lcc, kmax, j, k, i
   
   complex fre(wave_pts2)
   complex cc
!
   do i = 1, nlen
      cr(i) = u(i)
      cz(i) = 0.0
   end do
   call cfft(cr, cz, lnpt)
   call realtr(cr, cz, lnpt)
   do j = 1, nlen
      fre(j) = cmplx(cr(j), cz(j))
      cr(j) = 0.0
      cz(j) = 0.0
      u(j) = 0.0
   end do
!
! c1 = wave(pi2, 2), c2 = wave(pi, 2).
!
   u(1) = real(fre(2)*c1)
   u(2) = real(fre(2)*c2+fre(3)*c1)
   u(3) = real(-fre(2)*c2+fre(3)*c1)
   
   kmax = 2
   do j=3, jmax
      kmax = 2*kmax
      do k = 1, kmax
         cc = fre(k)*rwt1(k, j)+fre(k+kmax)*rwt2(k, j)
         cr(k) = real(cc)
         cz(k) = aimag(cc)
      end do
      lcc = j-1
      call cifft(cr, cz, Lcc)
      do k = 1, kmax
         u(kmax+k-1) = cr(k)
      end do
   end do
   end subroutine wavelet_obs


   pure subroutine wavelet_syn(cr, cz, u)
!
!  Args:
!  cr: real part of data vector
!  cz: imaginary part of data vector
!  u: wavelet coefficients from data vector
!
   implicit none
!       ******* Explaination ********
!       In order to normalized the amplitude, we move the coeffficient
!       sqrt(T/kmax)
!                         Jichen 1997, 1, 20
   real, intent(out) :: u(wave_pts2)
   real, intent(inout) :: cr(wave_pts2), cz(wave_pts2)
   integer lcc, kmax, lb, i, i1, i2, j, k
   real :: cr0(wave_pts2), cz0(wave_pts2)
!   complex fre(wave_pts2)
   complex cc

   lcc = nlen/2!2**(lnpt-1)
   lb = lcc-1
   lcc = lcc+1
   do i = 1, lb
      i1 = lcc+i
      i2 = lcc-i
      cr(i1) = cr(i2)
      cz(i1) = -cz(i2)
   end do
   cz(lcc) = 0.0
   do j = 1, nlen
      cr0(j) = cr(j)/nlen
      cz0(j) = cz(j)/nlen
   end do
!
! c1 = wave(pi2, 2), c2 = wave(pi, 2).
!
!   u(1) = real(fre(2)*c1)
!   u(2) = real(fre(2)*c2+fre(3)*c1)
!   u(3) = real(-fre(2)*c2+fre(3)*c1)
   u(1) = cr0(2)*real_c1 - cz0(2)*imag_c1
   u(2) = cr0(2)*real_c2 - cz0(2)*imag_c2
   u(2) = u(2) + cr0(3)*real_c1 - cz0(3)*imag_c1
   u(3) = -cr0(2)*real_c2 + cz0(2)*imag_c2
   u(3) = u(3) + cr0(3)*real_c1 - cz0(3)*imag_c1
        
   kmax = 2
   do j=3, jmax
      kmax = 2*kmax  
      do k = 1, kmax
         cr(k) = cr0(k)*wavelet_r1(k,j) - cz0(k)*wavelet_i1(k,j) &
         & + cr0(k+kmax)*wavelet_r2(k,j) - cz0(k+kmax)*wavelet_i2(k,j)
         cz(k) = cr0(k)*wavelet_i1(k,j) + cz0(k)*wavelet_r1(k,j) &
         & + cr0(k+kmax)*wavelet_i2(k,j) + cz0(k+kmax)*wavelet_r2(k,j)
         !cc = fre(k)*rwt1(k, j)+fre(k+kmax)*rwt2(k, j)
         !cr(k) = real(cc)
         !cz(k) = aimag(cc)
      end do
      lcc = j-1
      call cifft(cr, cz, Lcc)
      do k = 1, kmax
         u(kmax+k-1) = cr(k)
      end do
   end do

   end subroutine wavelet_syn


   pure subroutine cfft(xr, xi, n)
!
!  Args:
!  xr: real part of data vector
!  xi: imaginary part of data vector
!  n: length of data vector
!
!
! old version
!
   implicit none
   integer, intent(in) :: n
   real, intent(inout) :: xr(:), xi(:)
   integer k, i, ib, nb, lx, l, lb, lbh, ist, jh, j1, j
   real wkr, wki, qr, qi, holdr, holdi, flx
   LX = 2**N
   FLX = real(LX)
   NB = 1
   LB = LX
   DO L = 1, N
      LBH = LB/2
      DO IB = 1, NB           ! 2 ** (l-1) operaciones
         WKR = cos_fft(ib)
         WKI = -sin_fft(ib)
         IST = LB*(IB-1)
         DO J = IST+1, IST+LBH
            JH = J+LBH
            QR = XR(JH)*WKR-XI(JH)*WKI
            QI = XR(JH)*WKI+XI(JH)*WKR
            XR(JH) = (XR(J)-QR)
            XI(JH) = (XI(J)-QI)
            XR(J) = (XR(J)+QR)
            XI(J) = (XI(J)+QI)
         end do
      end do
      NB = 2*NB
      LB = LB / 2
   end do
   K = 0
   DO J = 1, LX
      K = KKK(J, N)
      IF(K.LT.J) cycle
      HOLDR = XR(J)
      HOLDI = XI(J)
      J1 = K+1
      XR(J) = XR(J1)
      XI(J) = XI(J1)
      XR(J1) = HOLDR
      XI(J1) = HOLDI
   end do
   DO I = 1, LX
      XR(I) = XR(I)/FLX
      XI(I) = XI(I)/FLX
   ENDDO
   END subroutine cfft


   pure subroutine cifft(xr, xi, n)
!
!  Args:
!  xr: real part of data vector
!  xi: imaginary part of data vector
!  n: length of data vector
!
   implicit none
!
! old version
!
   integer, intent(in) :: n
   real, intent(inout) :: xr(wave_pts2), xi(wave_pts2)
   integer k, ib, nb, lx, l, lb, lbh, ist, jh, j1, j
   real wkr, wki, qr, qi, holdr
   LX = 2**N
   NB = 1
   LB = LX
   DO L = 1, N
      LBH = LB/2
      DO IB = 1, NB           ! 2 ** n operations
         WKR = cos_fft(ib)
         WKI = sin_fft(ib)
         IST = LB*(IB-1)
         DO J = IST+1, IST+LBH
            JH = J+LBH
            QR = XR(JH)*WKR-XI(JH)*WKI
            QI = XR(JH)*WKI+XI(JH)*WKR
            XR(JH) = XR(J)-QR
            XI(JH) = XI(J)-QI
            XR(J) = XR(J)+QR
            XI(J) = XI(J)+QI
         ENDDO
      ENDDO
      NB = 2*NB
      LB = LB / 2
   ENDDO
   K = 0
   DO J = 1, LX
      K = KKK(j, N)
      IF(K.LT.J) cycle
      HOLDR = XR(J)
      J1 = K+1
      XR(J) = XR(J1)
      XR(J1) = HOLDR
   end do
   end subroutine cifft


   subroutine fourier_coefs()
   implicit none
!
! in this subroutine, we load into memory certain coeficients and values which are frequently used 
! in the method of the cfft. as the code spent a large amount of time in said part of cfft, we decided to
! load such data to memory to reduce computation time
!
   integer i, n, nb, nnb, ib, k
   real*8 :: omega

   do n = 2, 12
      k = 0
      nb = 2 ** (n-1)
      nnb = 2 ** n
      do ib = 1, nb
         omega = twopi*dble(k)/dble(nnb)
         cos_fft(ib) = cos(omega)
         sin_fft(ib) = sin(omega)
         do i = 2, n
            if (k .lt. 2**(n-i)) exit
            k = k-2**(n-i)
         end do
         k = k+2**(n-i)         
      end do
      k = 0
      do ib = 1, nnb
         kkk(ib, n) = k
         do i = 1, n
            if (k .lt. 2**(n-i)) exit
            k = k-2**(n-i)
         end do
         k = k+2**(n-i)         
      end do
   end do

   end subroutine fourier_coefs


   subroutine meyer_yamada()
   implicit none
   real*8 :: omega1, omega2
   integer j, is, kmax
   c1 = 2.*wave(twopi, 2)
   c2 = 2.*wave(pi, 2)
   real_c1 = real(c1)
   imag_c1 = aimag(c1)
   real_c2 = real(c2)
   imag_c2 = aimag(c2)
!
!       Create the coefficients of Mayer wavelet function
!       so it should be called before any further application.
!
   do j = 1, 12
      kmax = 2**(j-1)
      do is = 1, kmax
         omega1 = twopi*(is-1)/(kmax)
         omega2 = omega1+twopi
         rwt1(is, j) = 2.*wave(omega1, 2)
         rwt2(is, j) = 2.*wave(omega2, 2)
         wavelet_r1(is, j) = real(rwt1(is, j))
         wavelet_i1(is, j) = aimag(rwt1(is, j))
         wavelet_r2(is, j) = real(rwt2(is, j))
         wavelet_i2(is, j) = aimag(rwt2(is, j))
      end do
   end do
   end subroutine meyer_yamada


   function wave(w, ic) result(wave1)
!
!  Args:
!  w: frequency value
!  ic: 1 to get value of wavelet function, -1 if to get its complex conjugate
!
   implicit none
   real*8 :: w
   integer :: ic
   integer :: i, j, k
   real*8, parameter :: p43 = 4.d0*pi/3.d0, p23 = 2.d0*pi/3.d0
   real*8 :: g(2), p(2), fw, wf, wd, wp, wg, ww, hw, cct
   complex*16 :: pp, a
   complex :: wave1
   hw = w*0.5d0
   a = cmplx(cos(hw),-sin(hw), kind(1.d0))
   ww = abs(w)
   if (ww .gt. (2.d0*p43) .or. ww .lt. p23) then
      wave1 = cmplx(0.0, 0.0)
      return
   end if
   cct = 0.d0
   if ((w .ge. p23) .and. (w .le. p43)) then
      wd = w-p23
      wp = p43-w
      wg = exp(-1.d0/wd/wd)
      wf = exp(-1.d0/wp/wp)
      cct = exp(-0.5d0/wd/wd) / sqrt(wg + wf)
   elseif ((hw .ge. p23) .and. (hw .le. p43)) then
      wd = hw-p23
      wp = p43-hw
      wg = exp(-1.d0/wd/wd)
      wf = exp(-1.d0/wp/wp)
      cct = exp(-0.5d0/wp/wp) / sqrt(wg + wf)
   endif
   pp = a*cct
   if (ic .eq. 1) then
      wave1 = pp 
   else
      wave1 = cmplx(real(pp),-aimag(pp))    
   end if
   end function wave


   subroutine realtr(xr, xi, n)
!
!  Args:
!  xr: real part of data vector
!  xi: imaginary part of data vector
!  n: length of data vector
!
   implicit none
   real xr(:), xi(:)
   integer :: n
   integer :: i, i1, i2, lh, lb
   lh = 2**(n-1)
   lb = lh-1
   lh = lh+1
   do i = 1, lb
      i1 = lh+I
      i2 = lh-I
      xr(i1) = xr(i2)
      xi(i1) = -xi(i2)
   end do
   xi(lh) = 0.0
   end subroutine realtr


end module wavelets
