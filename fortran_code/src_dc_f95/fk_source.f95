!***************************************************************
!  F-K: @(#) source.f			1.0 4/29/2000
!
!  Copyright (c) 2000 by L. Zhu
!  See README file for copying and redistribution conditions.
!
! Setup the displacement-stress vector jump across the source
! The coefficients can be all derived from the potential jump in
! Haskell (1964), or directly from Takiuchi and Saito (1972),
! see ZR Appendix B.
!
!  Input:
!	stype	--- source type, 0=ex; 1=sf; 2=dc
!	xi	--- mu/(lambda+2*mu).
!	mu	--- shear modulus.
!  Output:
!	s(3,6)	--- source coef. for n=0, 1, 2
!
!  called by:	main()	in fk.f
!
!  modified history:
!	May 12, 2000	Lupei Zhu	initial coding.
!	July 17, 2000	Lupei Zhu	change xi, mu to be complex.
!	Oct. 30, 2008	Lupei Zhu	change all back to real.
!	Jan. 20, 2010   Lupei Zhu	add flip model up-side down.
!***************************************************************
module fk_source


   implicit none
   real*8 :: si(3, 6)


contains


   subroutine source(stype, xi, mu, flip)
!
!  Args:
!  stype: type of seismic source
!  xi: square of S velocity / square of P velocity at source
!  mu: shear modulous at source
!  flip: whether source is above or below receiver
!
   IMPLICIT NONE
   integer stype, flip
   real*8 :: xi, mu

   si(:,:) = 0.d0

   if (stype .EQ. 2) then                     ! DOUBLE_COUPLE
      si(1,2) = 2.d0*xi/mu
      si(1,4) = 4.d0*xi-3.d0
      si(2,1) = flip/mu
      si(2,5) =-si(2,1)
      si(3,4) = 1.d0
      si(3,6) =-1.d0
   else if (stype .EQ. 0) then                ! EXPLOSION
      si(1,2) = xi/mu
      si(1,4) = 2.d0*xi
   else if (stype .EQ. 1) then                ! SINGLE_FORCE, multiplied by k
      si(1,3) =-flip
      si(2,4) = -1.d0
      si(2,6) =  1.d0
   else                                       ! unknow source type
      write(0,*)'unknown source type'
      call exit(1)
   endif
!   write(0,*)flip, mu, xi
!   write(0,*)'source'
!   write(0,*)si

   end subroutine source


end module fk_source
