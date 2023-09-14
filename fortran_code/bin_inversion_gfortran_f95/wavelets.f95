module wavelets


   use constants, only : pi, twopi, wave_pts2
   use wavelet_param, only : get_data_param
   implicit none
   complex :: wavelet1(wave_pts2, 12), wavelet2(wave_pts2, 12), c1, c2
   integer :: switch_index(4200, 15)
   real :: cos_fft(4200), sin_fft(4200)
!   real :: wavelet_r1(wave_pts2, 12), wavelet_r2(wave_pts2, 12)
!   real :: wavelet_i1(wave_pts2, 12), wavelet_i2(wave_pts2, 12)
!   real :: real_c1, imag_c1, real_c2, imag_c2
   integer :: jmin, jmax, lnpt, nlen


contains
   

   subroutine wavelets_set_data_properties()
   implicit none
   integer :: max_freq
   call get_data_param(lnpt, jmin, jmax, nlen, max_freq)
   end subroutine wavelets_set_data_properties


   pure subroutine wavelet_syn(real1, imag1, coeffs)
!
!  Args:
!  real1: real part of data vector
!  imag1: imaginary part of data vector
!  coeffs: wavelet coefficients from data vector
!
   implicit none
!       ******* Explaination ********
!       In order to normalized the amplitude, we move the coeffficient
!       sqrt(T/kmax)
!                         Jichen 1997, 1, 20
   real, intent(out) :: coeffs(wave_pts2)
   real, intent(inout) :: real1(wave_pts2), imag1(wave_pts2)
   integer index1, n, kmax, half_length, length1, j1, j2, j, k
!   real :: real0(wave_pts2), imag0(wave_pts2)
   complex fft_array(wave_pts2)
   complex cc

   half_length = nlen/2!2**(lnpt-1)
   length1 = half_length-1
   index1 = half_length+1
   do j = 1, length1
      j1 = index1+j
      j2 = index1-j
      real1(j1) = real1(j2)
      imag1(j1) = -imag1(j2)
   end do
   imag1(index1) = 0.0
   do j = 1, nlen
      fft_array(j) = cmplx(real1(j), imag1(j))/nlen
!      real0(j) = real1(j)/nlen
!      imag0(j) = imag1(j)/nlen
!      real1(j) = 0.0
!      imag1(j) = 0.0
      coeffs(j) = 0.0
   end do
   
!   coeffs(1) = real0(2)*real_c1 - imag0(2)*imag_c1
!   coeffs(2) = real0(2)*real_c2 - imag0(2)*imag_c2
!   coeffs(2) = coeffs(2) + real0(3)*real_c1 - imag0(3)*imag_c1
!   coeffs(3) = -real0(2)*real_c2 + imag0(2)*imag_c2
!   coeffs(3) = coeffs(3) + real0(3)*real_c1 - imag0(3)*imag_c1
   coeffs(1) = real(fft_array(2)*c1)
   coeffs(2) = real(fft_array(2)*c2+fft_array(3)*c1)
   coeffs(3) = real(-fft_array(2)*c2+fft_array(3)*c1)
        
   kmax = 2
   do j=3, jmax
      kmax = 2*kmax  
      do k = 1, kmax
!         real1(k) = real0(k)*wavelet_r1(k,j) - imag0(k)*wavelet_i1(k,j) &
!         & + real0(k+kmax)*wavelet_r2(k,j) - imag0(k+kmax)*wavelet_i2(k,j)
!         imag1(k) = real0(k)*wavelet_i1(k,j) + imag0(k)*wavelet_r1(k,j) &
!         & + real0(k+kmax)*wavelet_i2(k,j) + imag0(k+kmax)*wavelet_r2(k,j)
         cc = fft_array(k)*wavelet1(k, j)+fft_array(k+kmax)*wavelet2(k, j)
         real1(k) = real(cc)
         imag1(k) = aimag(cc)
      end do
      n = j-1
      call cifft(real1, imag1, n)
      do k = 1, kmax
         coeffs(kmax+k-1) = real1(k)
      end do
   end do

   end subroutine wavelet_syn


   pure subroutine cfft(real1, imag1, n)
!
!  Args:
!  real1: real part of data vector
!  imag1: imaginary part of data vector
!  n: length of data vector
!
!
! old version
!
   implicit none
   integer, intent(in) :: n
   real, intent(inout) :: real1(:), imag1(:)
   integer k, i, chunks, length1, l, length2, length3, start, j1, j
   real real2, imag2, real3, imag3, norm
   length1 = 2**n
   norm = real(length1)
   chunks = 1
   length2 = length1
   DO L = 1, N
      length3 = length2/2
      DO i = 1, chunks           ! 2 ** (l-1) operaciones
         real2 = cos_fft(i)
         imag2 = -sin_fft(i)
         start = length2*(i-1)
         DO J = start+1, start+length3
            j1 = J+length3
            real3 = real1(j1)*real2-imag1(j1)*imag2
            imag3 = real1(j1)*imag2+imag1(j1)*real2
            real1(j1) = (real1(j)-real3)
            imag1(j1) = (imag1(j)-imag3)
            real1(j) = (real1(j)+real3)
            imag1(j) = (imag1(j)+imag3)
         end do
      end do
      chunks = 2*chunks
      length2 = length2 / 2
   end do
   K = 0
   DO J = 1, length1
      K = switch_index(J, N)
      IF(K.LT.J) cycle
      real2 = real1(J)
      imag2 = imag1(J)
      J1 = K+1
      real1(J) = real1(J1)
      imag1(J) = imag1(J1)
      real1(J1) = real2
      imag1(J1) = imag2
   end do
   DO I = 1, length1
      real1(I) = real1(I)/norm
      imag1(I) = imag1(I)/norm
   ENDDO
   END subroutine cfft


   pure subroutine cifft(real1, imag1, n)
!
!  Args:
!  real1: real part of data vector
!  imag1: imaginary part of data vector
!  n: length of data vector
!
   implicit none
!
! old version
!
   integer, intent(in) :: n
   real, intent(inout) :: real1(wave_pts2), imag1(wave_pts2)
   integer k, i, chunks, length1, l, length2, length3, start, j1, j
   real real2, imag2, real3, imag3
   length1 = 2**n
   chunks = 1
   length2 = length1
   DO L = 1, N
      length3 = length2/2
      DO i = 1, chunks           ! 2 ** n operations
         real2 = cos_fft(i)
         imag2 = sin_fft(i)
         start = length2*(i-1)
         DO J = start+1, start+length3
            j1 = J+length3
            real3 = real1(j1)*real2-imag1(j1)*imag2
            imag3 = real1(j1)*imag2+imag1(j1)*real2
            real1(j1) = real1(j)-real3
            imag1(j1) = imag1(j)-imag3
            real1(j) = real1(j)+real3
            imag1(j) = imag1(j)+imag3
         ENDDO
      ENDDO
      chunks = 2*chunks
      length2 = length2 / 2
   ENDDO
   K = 0
   DO J = 1, length1
      K = switch_index(j, N)
      IF(K.LT.J) cycle
      real2 = real1(J)
      J1 = K+1
      real1(J) = real1(J1)
      real1(J1) = real2
   end do
   end subroutine cifft


   subroutine fourier_coefs()
   implicit none
!
! in this subroutine, we load into memory certain coeficients and values which are frequently used 
! in the method of the cfft. as the code spent a large amount of time in said part of cfft, we decided to
! load such data to memory to reduce computation time
!
   integer i, j, n, length1, length2, k
   real*8 :: omega

   do n = 2, 12
      k = 0
      length1 = 2 ** (n-1)
      length2 = 2 ** n
      do i = 1, length1
         omega = twopi*dble(k)/dble(length2)
         cos_fft(i) = cos(omega)
         sin_fft(i) = sin(omega)
         do j = 2, n
            if (k .lt. 2**(n-j)) exit
            k = k-2**(n-j)
         end do
         k = k+2**(n-j)         
      end do
      k = 0
      do i = 1, length2
         switch_index(i, n) = k
         do j = 1, n
            if (k .lt. 2**(n-j)) exit
            k = k-2**(n-j)
         end do
         k = k+2**(n-j)         
      end do
   end do

   end subroutine fourier_coefs


   subroutine meyer_yamada()
   implicit none
   real*8 :: omega1, omega2
   integer j, i, kmax
   c1 = 2.*wave(twopi, 2)
   c2 = 2.*wave(pi, 2)
!   real_c1 = real(c1)
!   imag_c1 = aimag(c1)
!   real_c2 = real(c2)
!   imag_c2 = aimag(c2)
!
!       Create the coefficients of Mayer wavelet function
!       so it should be called before any further application.
!
   do j = 1, 12
      kmax = 2**(j-1)
      do i = 1, kmax
         omega1 = twopi*(i-1)/(kmax)
         omega2 = omega1+twopi
         wavelet1(i, j) = 2.*wave(omega1, 2)
         wavelet2(i, j) = 2.*wave(omega2, 2)
!         wavelet_r1(i, j) = real(wavelet1(i, j))
!         wavelet_i1(i, j) = aimag(wavelet1(i, j))
!         wavelet_r2(i, j) = real(wavelet2(i, j))
!         wavelet_i2(i, j) = aimag(wavelet2(i, j))
      end do
   end do
   end subroutine meyer_yamada


   function wave(omega, ic) result(wave1)
!
!  Args:
!  w: frequency value
!  ic: 1 to get value of wavelet function, else to get its complex conjugate
!
   implicit none
   real*8 :: omega
   integer :: ic
   real*8, parameter :: p43 = 4.d0*pi/3.d0, p23 = 2.d0*pi/3.d0
   real*8 :: fun1, fun2, abs_value, omega23, omega43, omega2, omega3
   complex*16 :: wave0, rotation
   complex :: wave1
   omega2 = omega*0.5d0
   rotation = cmplx(cos(omega2),-sin(omega2), kind(1.d0))
   omega3 = abs(omega)
   if (omega3 .gt. (2.d0*p43) .or. omega3 .lt. p23) then
      wave1 = cmplx(0.0, 0.0)
      return
   end if
   abs_value = 0.d0
   if ((omega .ge. p23) .and. (omega .le. p43)) then
      omega23 = omega-p23
      omega43 = p43-omega
      fun1 = exp(-1.d0/omega23/omega23)
      fun2 = exp(-1.d0/omega43/omega43)
      abs_value = exp(-0.5d0/omega23/omega23) / sqrt(fun1 + fun2)
   elseif ((omega2 .ge. p23) .and. (omega2 .le. p43)) then
      omega23 = omega2-p23
      omega43 = p43-omega2
      fun1 = exp(-1.d0/omega23/omega23)
      fun2 = exp(-1.d0/omega43/omega43)
      abs_value = exp(-0.5d0/omega43/omega43) / sqrt(fun1 + fun2)
   endif
   wave0 = rotation*abs_value
   if (ic .eq. 1) then
      wave1 = wave0 
   else
      wave1 = cmplx(real(wave0),-aimag(wave0))    
   end if
   end function wave


   subroutine realtr(real1, imag1, n)
!
!  Args:
!  real1: real part of data vector
!  imag1: imaginary part of data vector
!  n: length of data vector
!
   implicit none
   real real1(:), imag1(:)
   integer :: n
   integer :: i, i1, i2, length1, index1, half_length
   half_length = 2**(n-1)
   length1 = half_length-1
   index1 = half_length+1
   do i = 1, length1
      i1 = index1+I
      i2 = index1-I
      real1(i1) = real1(i2)
      imag1(i1) = -imag1(i2)
   end do
   imag1(index1) = 0.0
   end subroutine realtr


end module wavelets
