!***************************************************************
!  F-K: @(#) haskell.f			1.0 4/29/2000
!
!  Copyright (c) 2000 by L. Zhu
!  See README file for copying and redistribution conditions.
!
!  This file contains the following subroutines for computing Haskell propagator
!  matrix and other related matrices/vectors for the dynamic case.
!	layerParameter()	compute some parameters for each layer.
!       sh_ch()			compute Cosh(x), Sinh(x)/x, and x*Sinh(x).
!	haskell()		compute the P-SV Haskell matrix.
!       compoundMatrix()	compute the compound matrix of Haskell matrix.
!	eVector()		compute the compound matrix E|_12^ij.
!       initialG()		initialize g vector in the bottom half space.
!***************************************************************
module haskell


   use constants, only : one, two, nlay
   use vel_model_data, only : new_thick, mu
   use layer, only : kd, mu2, exa, exb, ra, rb, Ca, Cb, Ya, Xa, Yb, Xb, r, &
             & r1, y, y1, ka, kb
   implicit none


contains
   

   subroutine layerParameter(k, lay)!, kd, mu2, ra, rb, r, r1)
!***************************************************************
! compute some parameters for this layer
!	IN:
!		k   --- wave-number
!		lay --- layer number
!		volocity model passed in common/model/
!	Out:
!		common/layer/
! called by: kernel()		in kernel.f
!***************************************************************
   IMPLICIT NONE
   integer, intent(in) :: lay
   real*8, intent(in) :: k
   real*8 k2
   complex*16 kka, kkb
!   real*8, intent(out) :: kd, mu2
!   complex*16, intent(out) :: ra,rb,r,r1
!   complex*16, intent(in) :: ka(nlay),kb(nlay)
   k2 = k*k
   kka = ka(lay)/k2
   kkb = kb(lay)/k2
   r = two/kkb ! gamma
   kd = k*new_thick(lay)
   mu2 = two*mu(lay)
   ra = zsqrt(one - kka)
   rb = zsqrt(one - kkb)
   r1 = one - one/r
   end subroutine layerParameter


   subroutine sh_ch(c,y,x,ex,a,kd)
!
!  Args:
!  a: ?
!  kd: ?
!
!***************************************************************
! compute c=cosh(a*kd); y=sinh(a*kd)/a; x=sinh(a*kd)*a
! and multiply them by ex=exp(-Real(a*kd)) to supress overflow
!
! called by: compoundMatrix()		in compound.f
!***************************************************************
   IMPLICIT NONE
   complex*16, intent(out) :: c,y,x
   complex*16, intent(in) :: a
   real*8, intent(in) :: kd
   real*8 :: r,i
   real*8, intent(out) :: ex
   y = kd*a
   r = real(y)
   i = aimag(y)
   ex = exp(-r)
   y = 0.5d0*cmplx(cos(i),dsin(i),kind(1d0))
   x = ex*ex*conjg(y)
   c = y + x
   x = y - x
   y = x/a
   x = x*a
   end subroutine sh_ch


   subroutine haskellMatrix(a)!, exa, exb, mu2, Ca, Cb, Ya, Yb, Xa, Xb, r, r1)
!
!  Args:
!  a: 4x4 P-SV Haskell matrix for layer
!
!***************************************************************
! compute 4x4 P-SV Haskell a for the layer, the layer parameter
! is passed in by common /layer/.
!***************************************************************
   IMPLICIT NONE
   complex*16 a(5,5)
!   real*8 exa,exb,mu2
!   complex*16 Ca,Cb,Ya,Xa,Yb,Xb,r,r1

   Ca = Ca*exb
   Xa = Xa*exb
   Ya = Ya*exb
   Cb = Cb*exa
   Yb = Yb*exa
   Xb = Xb*exa

! p-sv, scaled by exa*exb, see p381/Haskell1964 or EQ 17 of ZR
   a(1,1) = r*(Ca-r1*Cb)
   a(1,2) = r*(r1*Ya-Xb)
   a(1,3) = (Cb-Ca)*r/mu2
   a(1,4) = (Xb-Ya)*r/mu2

   a(2,1) = r*(r1*Yb-Xa)
   a(2,2) = r*(Cb-r1*Ca)
   a(2,3) = (Xa-Yb)*r/mu2
   a(2,4) =-a(1,3)

   a(3,1) = mu2*r*r1*(Ca-Cb)
   a(3,2) = mu2*r*(r1*r1*Ya-Xb)
   a(3,3) = a(2,2)
   a(3,4) =-a(1,2)

   a(4,1) = mu2*r*(r1*r1*Yb-Xa)
   a(4,2) =-a(3,1)
   a(4,3) =-a(2,1)
   a(4,4) = a(1,1)

! sh, the Haskell matrix is not needed. it is replaced by exb
   a(5,5) = exb

   end subroutine haskellMatrix


   subroutine compoundMatrix(a)!, exa, exb, Ca, Cb, Ya, Yb, Xa, Xb, kd, mu2, ra, rb, r, r1)
!
!  Args:
!  a: compound matrix of the P-SV Haskell matrix
!
!***************************************************************
! The upper-left 5x5 is the 6x6 compound matrix of the P-SV Haskell matrix,
!	a(ij,kl) = A|_kl^ij, ij=12,13,14,23,24,34,
! after dropping the 3rd row and colummn and replacing the 4th row
! by (2A41, 2A42, 2A44-1,2A45,2A46) (see W&H, P1035).
! The lower-right c 2x2 is the SH part of the Haskell matrix.
! Input: layer parameters passed in by /layer/.
! Output: compound matrix a, scaled by exa*exb for the P-SV and exb for the SH.
!***************************************************************
   IMPLICIT NONE
   real*8 ex
   complex*16, intent(out) :: a(7,7)
   complex*16 CaCb, CaXb, CaYb, XaCb, XaXb, YaCb, YaYb, r2, r3
!   real*8, intent(in) :: kd,mu2
!   real*8, intent(out) :: exa,exb
!   complex*16, intent(in) :: ra,rb,r,r1
!   complex*16, intent(out) :: Ca,Cb,Ya,Xa,Yb,Xb

   call sh_ch(Ca,Ya,Xa,exa,ra,kd)
   call sh_ch(Cb,Yb,Xb,exb,rb,kd)
   CaCb=Ca*Cb
   CaYb=Ca*Yb
   CaXb=Ca*Xb
   XaCb=Xa*Cb
   XaXb=Xa*Xb
   YaCb=Ya*Cb
   YaYb=Ya*Yb
   ex = exa*exb
   r2 = r*r
   r3 = r1*r1

! p-sv, scaled by exa*exb to supress overflow
!
   a(1,1) = ((one+r3)*CaCb-XaXb-r3*YaYb-two*r1*ex)*r2
   a(1,2) = (XaCb-CaYb)*r/mu2
   a(1,3) = ((one+r1)*(CaCb-ex)-XaXb-r1*YaYb)*r2/mu2
   a(1,4) = (YaCb-CaXb)*r/mu2
   a(1,5) = (two*(CaCb-ex)-XaXb-YaYb)*r2/(mu2*mu2)

   a(2,1) = (r3*YaCb-CaXb)*r*mu2
   a(2,2) = CaCb
   a(2,3) = (r1*YaCb-CaXb)*r
   a(2,4) =-Ya*Xb
   a(2,5) = a(1,4)

   a(3,1) = two*mu2*r2*(r1*r3*YaYb-(CaCb-ex)*(r3+r1)+XaXb)
   a(3,2) = two*r*(r1*CaYb-XaCb)
   a(3,3) = two*(CaCb - a(1,1)) + ex
   a(3,4) =-two*a(2,3)
   a(3,5) =-two*a(1,3)

   a(4,1) = mu2*r*(XaCb-r3*CaYb)
   a(4,2) =-Xa*Yb
   a(4,3) =-a(3,2)/two
   a(4,4) = a(2,2)
   a(4,5) = a(1,2)

   a(5,1) = mu2*mu2*r2*(two*(CaCb-ex)*r3-XaXb-r3*r3*YaYb)
   a(5,2) = a(4,1)
   a(5,3) =-a(3,1)/two
   a(5,4) = a(2,1)
   a(5,5) = a(1,1)

! sh, scaled by exb
   a(6,6) = Cb
   a(6,7) =-two*Yb/mu2
   a(7,6) =-mu2*Xb/two
   a(7,7) = Cb

   end subroutine compoundMatrix


   subroutine eVector(e)!, ra, rb, r1, mu2)
!
!  Args:
!  e: comes from E the similarity matrix 
!
!***************************************************************
! The first 5 members are E|_12^ij, ij=12,13,23,24,34.
! The last two are the first column of SH E matrix.
!***************************************************************
   IMPLICIT NONE
!   real*8 mu2
!   complex*16 ra, rb, r1
   complex*16 e(7)
! For p-sv, compute E|_(12)^(ij), ij=12, 13, 23, 24, 34.
   e(1) = ra*rb-one
   e(2) = mu2*rb*(one-r1)
   e(3) = mu2*(r1-ra*rb)
   e(4) = mu2*ra*(r1-one)
   e(5) = mu2*mu2*(ra*rb-r1*r1)
! sh part
   e(6) =-one
   e(7) = mu2*rb/two
   end subroutine eVector


   subroutine initialG(g)!, ra, rb, r, r1, mu2)
!
!  Args:
!  g: vector to be propagated
!
!***************************************************************
! Initialize the g row-vector. The first 5 elements are the
! inverse(E)|_{ij}^{12}, ij=12,13,23,24,34.
! g14 is omitted because g14 = -g23.
! The last two are the 5th row of E^-1.
!***************************************************************
   IMPLICIT NONE
   complex*16 g(7), delta
!   real*8 mu2
!   complex*16 ra, rb, r, r1
! p-sv, see EQ 33 on ZR/p623, constant omitted.
   delta  = r*(one-ra*rb)-one
   g(1) = mu2*(delta-r1)
   g(2) = ra
   g(3) = delta
   g(4) =-rb
   g(5) = (one+delta)/mu2
! sh, use the 5th row of E^-1, see EQ A4 on ZR/p625, 1/2 omitted
   g(6) =-one
   g(7) = two/(rb*mu2)
   end subroutine initialG


end module haskell
