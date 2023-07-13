!***************************************************************
!  F-K: @(#) kernel.f			1.0 4/29/2000
!
!  Copyright (c) 2000 by L. Zhu
!  See README file for copying and redistribution conditions.
!
! Compute dynamic displacement kernels from a point source in a multi-layer medium
! using the Haskell propagator matrix.
!
!	reference:
!		Haskell (1964), BSSA
!		Wang and Herrmann (1980), BSSA
!
!	with some modifications:
!	(1) modify the definition of B and K and put back w^2 so all matrix
!	    elements are either dimensionless or in in shear modular unit.
!	(2) suppress overflow by factoring out e^(Re(ra+rb)*d) in matrix
!	    elements.
!	(3) the cordinate is Z up, R outward and T counterclockwise
! Input:
!	k	wavenumber
!	/model/	velocity model (see model.h)
! Output:
!	u(i,3)	kernels for azimuthal mode i = 0, 1, 2 in the order of Z, R, T.
!		Note that 1/2*pi*mu is omitted.
! Called by
!	main() in fk.f.
!
! Subroutines called are in haskell.f (st_haskell.f in static case) and prop.f.
!***************************************************************
module fk_kernel


   use constants, only : nlay, zero, one, two, epsilon
   use vel_model_data, only : n_layers_new, src, rcv, mu, xi, new_thick
   use fk_source, only : si
   use layer, only : kd, mu2, exa, exb, ra, rb, Ca, Cb, Ya, Xa, Yb, Xb, r, r1, y, y1
   use haskell, only : layerParameter, haskellMatrix, compoundMatrix, eVector, initialG
   use prop, only : propagateG, initialZ, propagateZ, initialB, propagateB
   implicit none


contains


   subroutine kernel(k, u)!, ka, kb)
!
!  Args:
!  k: wavenumber
!  u: kernel of GF to be computed
!
   IMPLICIT NONE
!   real*8 kd,exa,exb,mu2
   integer stype,updn
!   complex*16 ra, rb, r, r1, Ca, Cb, Xa, Xb, Ya, Yb
   integer :: i, j
   real*8 :: k
!   complex*16 :: ka(nlay), kb(nlay)
!   complex u(3,3)
   complex*16 u(3,3)
   complex*16 rayl, love, dum
   complex*16 a(5,5), b(7,7), c(7,7), e(7), g(7), z(3,5), ss(3,6)
   stype = 2
   updn = 0
! Explanations:
! a --- 4x4 p-sv Haskell matrix and a(5,5)=exb.
! b --- product of compound matrices from the receiver to the surface.
! c --- compound matrix.
! e --- vector, the first 5 members are E0|_{12}^{ij}, ij=12,13,23,24,34;
!	the last two are the 1st column of the 2x2 SH E0 (or unit vector if the top is elastic).
! g --- vector containing the Rayleigh and Love denominators. It is initialized in the
!	bottom half-space with  (E^-1)|_{ij}^{12}, ij=12,13,23,24,34, and
!	the 1st row of the 2x2 SH E^-1 (or a unit vector if the bottom is vacume).
! z --- z(n,j)=s(i)*X|_{ij}^{12} for p-sv and s(i)*X_5i for sh.

   call initialB(b, e, g)
! propagation - start from the bottom
   do j = n_layers_new, 1, -1
      call layerParameter(k, j)!, ka, kb)!, kd, mu2, ra, rb, r, r1)
      if ( j.EQ.n_layers_new .AND. new_thick(j).LT.epsilon ) then ! half space
         call initialG(g)!, ra, rb, r, r1, mu2)
      else if ( j.EQ.1 .AND. new_thick(1).LT.epsilon ) then 
         call eVector(e)!, ra, rb, r1, mu2)
         exit
      else
         call compoundMatrix(c)!, exa, exb, Ca, Cb, Ya, Yb, Xa, Xb, kd, mu2, ra, rb, r, r1)
         call propagateG(c, g)
      endif
      if ( j.EQ.src ) then
         call separatS(ss)!, mu2, ra, rb, r, r1)
         call initialZ(ss, g, z)
      endif
      if ( j.LT.src ) then
         if ( j.GE.rcv ) then
            call haskellMatrix(a)!, exa, exb, mu2, Ca, Cb, Ya, Yb, Xa, Xb, r, r1) 
            call propagateZ(a, z)
         else
            call propagateB(c, b)
         endif
      endif
   enddo

! add the top halfspace boundary condition
   e(3) = two*e(3)
   rayl = g(1)*e(1)+g(2)*e(2)+g(3)*e(3)+g(4)*e(4)+g(5)*e(5)
   love = g(6)*e(6)+g(7)*e(7)
   do i=1, 4
      g(i) = zero
      do j=1, 5
         g(i) = g(i) + b(i,j)*e(j)
      enddo
   enddo
   g(3) = g(3)/two
   g(6) = b(6,6)*e(6) + b(6,7)*e(7) 
   do i=1, 3
      dum    = z(i,2)*g(1)+z(i,3)*g(2)-z(i,4)*g(3)
      z(i,2) =-z(i,1)*g(1)+z(i,3)*g(3)+z(i,4)*g(4)
      z(i,1) = dum
      z(i,5) = z(i,5)*g(6)
   enddo

! displacement kernels at the receiver
   dum = k
   if ( stype.EQ.1 ) dum = one
   do i = 1, 3
      u(1,i) = dum*z(i,2)/rayl
      u(2,i) = dum*z(i,1)/rayl
      u(3,i) = dum*z(i,5)/love
   enddo

   end subroutine kernel


   subroutine separatS(ss)!, mu2, ra, rb, r, r1)
!
!  Args:
!  ss: source vector (?) 
!
   IMPLICIT NONE
!   complex*16 ra, rb, r, r1
!   real*8 mu2
   integer stype,updn
   integer i, j, ii, jj
   complex*16 temp(4,4), temp_sh, ss(3,6), ra1, rb1, dum

   stype = 2
   updn = 0
   ii = stype+1
   if (updn.EQ.0) then
      ss = si
      return
   endif
!
! down-going (updn=1) or up-going(updn=-1) matrix: E*diag(...)*inv(E), without the 1/2 factor
!
   ra1 = one/ra
   rb1 = one/rb
   dum = updn*r
   temp(1,1) = one
   temp(1,2) = dum*(rb-r1*ra1)
   temp(1,3) = zero
   temp(1,4) = dum*(ra1-rb)/mu2
   temp(2,1) = dum*(ra-r1*rb1)
   temp(2,2) = one
   temp(2,3) = dum*(rb1-ra)/mu2
   temp(2,4) = zero
   temp(3,1) = zero
   temp(3,2) = dum*(rb-r1*r1*ra1)*mu2
   temp(3,3) = one
   temp(3,4) = dum*(r1*ra1-rb)
   temp(4,1) = dum*(ra-r1*r1*rb1)*mu2
   temp(4,2) = zero
   temp(4,3) = dum*(r1*rb1-ra)
   temp(4,4) = one
   temp_sh   = (updn*two/mu2)*rb1
!
   do i=1, ii
      do j=1,4
         dum = zero
         do jj = 1,4
            dum = dum + temp(j,jj)*si(i,jj)
         enddo
         ss(i,j) = dum/two
      enddo
      ss(i,5) = (si(i,5) + temp_sh*si(i,6))/two
      ss(i,6) = (si(i,6) + si(i,5)/temp_sh)/two
   enddo
   end subroutine separatS


end module fk_kernel
