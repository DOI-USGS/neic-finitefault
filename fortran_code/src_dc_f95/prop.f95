!***************************************************************
!  F-K: @(#) prop.f			1.0 4/29/2000
!
!  Copyright (c) 2000 by L. Zhu
!  See README file for copying and redistribution conditions.
!
!  This file contains the following subroutines for propagating Haskell matrices.
!       propagateG()		propoagate the g vector by the compound matrix.
!       initialZ()		initialize z vector at the source.
!	propagateZ()		propagate the z vector by mulitiplying the Haskell matrix.
!	initialB()		initialize the B matrix and  e/g vectors.
!	propagateB()		compound matrix product.
!***************************************************************


module prop


   use constants, only : zero, one
   implicit none


contains


   subroutine propagateG(a, g)
!
!  Args:
!  a: compound Haskell matrix
!  g: vector to be propagated
!
!***************************************************************
! propagate g vector upward using the compound matrix
!	g = g*a
!***************************************************************
   IMPLICIT NONE
!   include 'constants.h'
   integer :: i, j
   complex*16 :: g(7), a(7,7), temp(7)
! p-sv
   do i = 1, 5
      temp(i) = zero
      do j = 1, 5
         temp(i) = temp(i) + g(j)*a(j,i)
      enddo
   enddo
   do i = 1, 5
      g(i) = temp(i)
   enddo
! sh
   do i = 6, 7
      temp(i) = g(6)*a(6,i) + g(7)*a(7,i)
   enddo
   do i = 6, 7
      g(i) = temp(i)
   enddo
   end subroutine propagateG


   subroutine initialZ(s, g, z)
!
!  Args:
!  s: source vector
!  g: vector to be propagated
!  z: propagator matrix from the surface to above the source
!
!***************************************************************
! initialize the row-vector z at the source z(j)=s(i)*X|_ij^12
! for P-SV and z(j)=s(i)*X(5,i) for SH.
!  input:
!	s(3,6)	---- source coef. for n=0,1,2
!	g(7)	---- g vector used to construct matrix X|_ij^12
!		     |	0   g1  g2 -g3 |
!	 X|_ij^12 =  | -g1  0   g3  g4 | for P-SV.
!		     | -g2 -g3  0   g5 |
!		     |  g3 -g4 -g5  0  |
!	 X(5,i) = ( g6 g7 )	for SH.
!  output:
!	z(3,5)  ---- z vector for n=0,1,2
!***************************************************************
   IMPLICIT NONE
   integer :: i
   complex*16 :: g(7), z(3,5), s(3,6)
   do i=1, 3
! for p-sv, see WH p1018
      z(i,1) =-s(i,2)*g(1)-s(i,3)*g(2)+s(i,4)*g(3)
      z(i,2) = s(i,1)*g(1)-s(i,3)*g(3)-s(i,4)*g(4)
      z(i,3) = s(i,1)*g(2)+s(i,2)*g(3)-s(i,4)*g(5)
      z(i,4) =-s(i,1)*g(3)+s(i,2)*g(4)+s(i,3)*g(5)
! for sh
      z(i,5) = s(i,5)*g(6)+s(i,6)*g(7)
   enddo
   end subroutine initialZ


   subroutine propagateZ(a, z)
!
!  Args:
!  a: 4x4 P-SV Haskell matrix for layer
!  z: propagator matrix from the surface to above the source
!
!***************************************************************
!  apply the Haskell matrix a to the z vector
!	z = z*a
!***************************************************************
   IMPLICIT NONE
!   include 'constants.h'
   integer :: i, j, l
   complex*16 z(3,5), a(5,5), temp(4)
   do i = 1, 3
! p-sv
      do j = 1, 4
         temp(j) = zero
         do l = 1, 4
            temp(j) = temp(j) + z(i,l)*a(l,j)
         enddo
      enddo
      do j = 1, 4
         z(i,j) = temp(j)
      enddo
! sh, due to the exb scaling of the sh haskell matrix.
      z(i,5) = z(i,5)*a(5,5)
   enddo
   end subroutine propagateZ


   subroutine initialB(b, e, g)
!
!  Args:
!  b: 7x7 identity matrix
!  g: vector to be propagated
!  e: comes from E the similarity matrix 
!
!***************************************************************
! Initialize b as an unit matrix; e as an unit vector;
! e = (1 0 0 0 0 1 0) for top halfspace boundary condition;
! g = (0 0 0 0 1 0 1) for free bottom boundary condition.
!***************************************************************
   IMPLICIT NONE
!   include 'constants.h'
   integer :: i, j
   complex*16 b(7,7), e(7), g(7)
   do i=1, 7
      do j=1, 7
         b(i,j) = zero
      enddo
      b(i,i) = one
      e(i) = zero
      g(i) = zero
   enddo
   e(1) = one
   g(5) = one
   e(6) = one
   g(7) = one
   end subroutine initialB


   subroutine propagateB(c, b)
!
!  Args:
!  b: 7x7 identity matrix
!  c: compound matrix
!
!***************************************************************
!	b = b*c
!***************************************************************
   IMPLICIT NONE
!   include 'constants.h'
   integer i, j, l
   complex*16 b(7,7), c(7,7), temp(7)
! p-sv
   do i = 1, 5
      do j = 1, 5
         temp(j) = zero
         do l = 1, 5
            temp(j) = temp(j) + b(i,l)*c(l,j)
         enddo
      enddo
      do j = 1,5
         b(i,j) = temp(j)
      enddo
   enddo
! sh
   do i = 6, 7
      do j = 6, 7
         temp(j) = b(i,6)*c(6,j) + b(i,7)*c(7,j)
      enddo
      do j = 6, 7
         b(i,j) = temp(j)
      enddo
   enddo
   end subroutine propagateB


end module prop
