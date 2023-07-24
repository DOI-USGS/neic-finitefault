module random_gen


   implicit none
   integer, parameter, private :: ia=16807, im=2147483647, iq=127773, ir=2836, ntab=32, ndiv=1+int((im-1)/ntab)
   real, parameter, private :: am=1./im, eps=1.2e-7, rnmx=1.0-eps
   integer :: iv(ntab), iy
!   data iv/ntab*0/, iy/0/
   integer, parameter, private :: mbig=1000000000, midum=161803398, mz=0
   real, parameter, private :: fac=1.e-9
   integer :: iff, inext, inextp, ma(55)
!   data iff/0/
   integer :: seed


contains


   subroutine start_seed(idum)
!
!  Args:
!  idum: seed value
!
   implicit none
   integer, intent(in) :: idum
   seed = idum
   iff = 0
   iv(:) = 0
   iy = 0
   end subroutine start_seed

   
   function ran1() result(random)
   implicit none
   integer :: j, k
   real :: random
   if (seed .le. 0 .or. iy .eq. 0) then
      seed = max(-seed, 1)
      do j = ntab+8, 1, -1
         k = seed/iq
         seed = ia*(seed-k*iq)-ir*k
         if (seed .lt. 0) seed = seed+im
         if (j .le. ntab) iv(j) = seed
      end do
      iy = iv(1)
   end if
   k = seed/iq
   seed = ia*(seed-k*iq)-ir*k
   if (seed .lt. 0) seed = seed+im
   j = 1+iy/ndiv
   iy = iv(j)
   iv(j) = seed
   random = min(am*iy, rnmx)
   end function ran1


   function ran3() result(random)
!           
!*****  routine to generate a uniformly distributed random *****
!*****  number on the interval [0, 1].                      *****
!
   implicit none
   integer :: i, ii, k, mj, mk
   real :: random
   if (seed .lt. 0 .or. iff .eq. 0) then
      iff = 1
      mj = midum-iabs(seed)
      mj = mod(mj, mbig)
      ma(55) = mj
      mk = 1
      do i = 1, 54
         ii = mod(21*i, 55)
         ma(ii) = mk
         mk = mj-mk
         if (mk .lt. mz) mk = mk+mbig
         mj = ma(ii)
      end do
      do k = 1, 4
         do i = 1, 55
            ma(i) = ma(i)-ma(1+mod(i+30, 55))
            if (ma(i) .lt. mz) ma(i) = ma(i)+mbig
         end do
      end do
      inext = 0
      inextp = 31
      seed = 1
   end if
   inext = inext+1
   if (inext .eq. 56) inext = 1
   inextp = inextp+1
   if (inextp .eq. 56) inextp = 1
   mj = ma(inext)-ma(inextp)
   if (mj .lt. mz) mj = mj+mbig
   ma(inext) = mj
   random = mj*fac
   end function ran3


   subroutine cauchy(t, x)
!
!  Args:
!  i: temṕerture of annealing method
!  x: output
!
   implicit none
   real, intent(in) :: t
   real, intent(out) :: x
   real :: u, sgn, uu
   u = ran3()
   sgn = 1.
   if ((u-0.5) .lt. 0.) sgn = -1.
   uu = abs(2.*u-1.)
   x = sgn*t*((1.+1./t)**uu-1.)
   end subroutine cauchy


end module random_gen
