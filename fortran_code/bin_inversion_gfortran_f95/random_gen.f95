module random_gen
!
! TODO: there is an error in ran3. This routine can, sometimes, output negative numbers!
! TODO: one case where this happens, is when the input seed value for ran3 is 1731825384
! TODO: why not replace ran1, ran3 with gfortran's random_number? 
!


   implicit none
   integer, parameter, private :: constant1=16807, big_constant1=2147483647, constant3=127773
   integer, parameter, private :: constant4=2836, size1=32, ndiv=1+int((big_constant1-1)/size1)
   real, parameter, private :: factor1=1./big_constant1, eps=1.2e-7, rnmx=1.0-eps
   integer :: random_list1(size1), random_int
!   data random_list1/size1*0/, random_int/0/
   integer, parameter, private :: big_constant2=1000000000, mseed=161803398
   real, parameter, private :: factor2=1.e-9
   integer :: initialize, index1, index2, random_list2(55)
!   data initconstant1lize1/0/
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
   initialize = 0
   random_list1(:) = 0
   random_int = 0
   end subroutine start_seed

   
   function ran1() result(random)
   implicit none
   integer :: j, k
   real :: random
   if (seed .le. 0 .or. random_int .eq. 0) then
      seed = max(-seed, 1)
      do j = size1+8, 1, -1
         k = seed/constant3
         seed = constant1*(seed-k*constant3)-constant4*k
         seed = modulo(seed, big_constant1)
!         if (seed .lt. 0) seed = seed+big_constant1
         if (j .le. size1) random_list1(j) = seed
      end do
      random_int = random_list1(1)
   end if
   k = seed/constant3
   seed = constant1*(seed-k*constant3)-constant4*k
   seed = modulo(seed, big_constant1)
!   if (seed .lt. 0) seed = seed+big_constant1
   j = 1+random_int/ndiv
   random_int = random_list1(j)
   random_list1(j) = seed
   random = min(factor1*random_int, rnmx)
!   if (random .lt. 0) write(*,*)'ran1 LT 0: ', random
!   if (random .gt. 1) write(*,*)'ran1 GT 1: ', random
   end function ran1


   function ran3() result(random)
!           
!*****  routine to generate a uniformly distributed random *****
!*****  number on the interval [0, 1].                      *****
!
!*****************************************************************************
!   This yiels an erroneous negative value when the input seed is 1731825384! 
!*****************************************************************************
!
   implicit none
   integer :: i, j, k, int1, int2
   real :: random
   if (seed .lt. 0 .or. initialize .eq. 0) then
      initialize = 1
      int1 = mseed-abs(seed)
      int1 = mod(int1, big_constant2)
!      write(*,*)seed, int1
      random_list2(55) = int1
      int2 = 1
      do i = 1, 54
         j = mod(21*i, 55)
         random_list2(j) = int2
         int2 = int1-int2
!         int2 = modulo(int2, big_constant2)
         if (int2 .lt. 0) int2 = int2+big_constant2
         int1 = random_list2(j)
      end do
      do k = 1, 4
         do i = 1, 55
            random_list2(i) = random_list2(i)-random_list2(1+mod(i+30, 55))
!            random_list2(i) = modulo(random_list2(i), big_constant2) 
            if (random_list2(i) .lt. 0) random_list2(i) = random_list2(i)+big_constant2
         end do
      end do
!      write(*,*)seed
!      write(*,*)random_list2
      index1 = 0
      index2 = 31
      seed = 1
   end if
   index1 = modulo(index1, 55)
   index1 = index1+1
!   if (index1 .eq. 56) index1 = 1
   index2 = modulo(index2, 55)
   index2 = index2+1
!   if (index2 .eq. 56) index2 = 1
   int1 = random_list2(index1)-random_list2(index2)
!   int1 = modulo(int1, big_constant2)
   if (int1 .lt. 0) int1 = int1+big_constant2
   random_list2(index1) = int1
   random = int1*factor2
!   if (random .lt. 0) write(*,*)'ran3 LT 0: ', random
!   if (random .gt. 1) write(*,*)'ran3 GT 1: ', random
   end function ran3


   subroutine cauchy(temp, random)
!
!  Args:
!  temp: temperature of annealing method
!  random: output
!
   implicit none
   real, intent(in) :: temp
   real, intent(out) :: random
   real :: value1, sign1, value2
   !call random_number(value1)
   value1 = ran3()
   !sign1 = 1.
   !if (value1 .lt. 0.5) sign1 = -1.
   sign1 = sign(1., value1 - 0.5)
   value2 = abs(2.*value1-1.)
   random = sign1*temp*((1.+1./temp)**value2-1.)
   end subroutine cauchy


end module random_gen
