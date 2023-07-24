!**********************************************************************
!  F-K:	@(#) fk.f			3.0 01/20/2011
!
!  Copyright (c) 1999-2011 by L. Zhu
!  See README file for copying and redistribution conditions.
!
! Compute displacements from a point source in layered
! media using the Thompson-Haskell propagator matrix.
! The source can have azimuthal mode number up to n=2, including explosion,
! single force, and double-couple. All source time functions are Dirac delta.
!
! The output displacements are in SAC format in the order of
! vertical (UP), radial, and tangential (couterclockwise) for n=0, 1, 2
! (i.e. Z0, R0, T0, Z1, ...).
! Their units are (assumme v in km/s, rho in g/cm^3, thickness in km):
!	10^-20 cm/(dyne cm)	for double couple source and explosion
!	10^-15 cm/dyne		for single force
!
! The code uses summmation for wave-number (k) integration and then FFT for
! omega integration.
!
! subroutines called:
!	kernel() in kernel.f	computing Kernels U(k,w)
!
! Written by Lupei Zhu, Seismo Lab, Caltech (lupei@gps.caltech.edu)
!
! Modified history:
!	03/05/1996  Lupei Zhu	Initial coding.
!	03/07/1999  Lupei Zhu	Use conditional compiling option
!				to consolidate the options of
!				exact or approx. Bessel functions.
!	05/06/1999  Lupei Zhu	Remove the division of (i*w) which
!				was put for getting the step-response
!				of the double-couple source.
!	04/29/2000  Lupei Zhu   Determine kmax using source depth hs.
!	05/06/2000  Lupei Zhu	Use total and count to estimate progress.
!	05/06/2000  Lupei Zhu	Use array instead of writing temp to disk.
!	07/17/2000  Lupei Zhu	Treat mu as complex and move source()
!				back into the w-k loops.
!	07/26/2000  Lupei Zhu	Modify to do sigma correction on real
!				array. The previous one does this on
!				complex array and causes long-period noise.
!	05/02/2005  Lupei Zhu	Multiply i*w to get impulse resp. for the
!				single-force source.
!	07/19/2005  Lupei Zhu   Add a hi-pass filter (wc1,wc2).
!	11/07/2006  Lupei Zhu	Offset comp. indices of outputs for explosion.
!	10/29/2008  Lupei Zhu	Add non-free surface boundary condition and
!				receiver interface, also change mu to real
!				and move source() out of the w-k loops.
!	11/06/2008  Lupei Zhu	Combine the static and dynamic cases.
!	10/05/2010  Lupei Zhu	Add option for up- or down-going wave (updn)
!	01/20/2011  Lupei Zhu	Correct two bugs (flip and phase shift)
!				for single force.
!				Correct a bug in free bottom boundary condi.
!				Use thickness<epsilon to indicate halfspace.

!                   pk          Added small modifications.
!**********************************************************************
!      program fk
module fk_openmp


   use constants, only : nlay, twopi, pi, ndis, nt
   use vel_model_data, only : depths, xi, mu, new_vel_p, new_vel_s, rcv, src, new_dens, qqp, qqs, flip_model, &
                        &     extra, extra2, n_layers_new
   use layer, only : ka, kb
   use fk_source, only : source
   use fk_kernel, only : kernel
   use retrieve_gf, only : lnpt, dt, d_step
   use bessel2, only : aj0s, aj1s, aj2s
   implicit none
   integer, parameter, private :: double=kind(1d0)


contains


   subroutine sub_bs_dc(nx,x,t0,green,disp)
!
!  Args:
!  nx: amount of source-station distances to compute GF
!  x: range of source-station distances to compute GF
!  t0: start time for each waveform
!  green: GF to be computed
!  disp: True if GF for cGPS, False otherwise
!
!
   use omp_lib
   IMPLICIT NONE
   integer stype,updn
   complex*16 :: att,w,w2
   logical dynamic, disp
!
   real :: green(nt, 8, ndis)
!
   integer i,j,l,nfft,nfft2,nfft3,n,ix,nx,tenPerc,count,total
   integer nCom,wc,wc1,wc2,smth,idx0,flip,new_old(9)
   real*8 :: k,omega,dk,dw,sigma,const,phi,hs,xmax,vs
   real*8 :: dfac,pmin,pmax,kc,taper,filter,z2
   real x(ndis),t0(ndis)
!
   integer ixx
!
   complex*16 :: nf1,nf2,u(3,3), aux(9, ndis)
!   complex :: nf1,nf2,u(3,3), aux(9, ndis)
   real aj0,aj1,aj2,z,tdata(2*nt)
!
   complex, allocatable :: summ(:, :, :)
   complex data1(nt), kahan_y(9, ndis), kahan_t(9, ndis), kahan_c(9, ndis)
   allocate(summ(9,ndis,nt))

   call omp_set_num_threads(48)
   write(0,*)'Compute near field GF at fixed depth...'
   dynamic = .TRUE.
   nCom = 9
!       sequence 1   2   3  4    5   6   7   8
!               tss tds rss rds r45 zss zds z45
   new_old(1)=9
   new_old(2)=6
   new_old(3)=8
   new_old(4)=5
   new_old(5)=2
   new_old(6)=7
   new_old(7)=4
   new_old(8)=1
!   call cmodel_fk(depth,0.0, jjo, nlay_s, nlay_r, a, b,
!     *  den, th, qa, qb)

   stype=2
   sigma=2.d0
   nfft=2**lnpt
   taper=0.2d0
   smth=1
   updn=0
   wc1=1
   wc2=1

!***************************************************************
! input velocity model
!      write(0,'(a)') 'Input nlay src_layer src_type rcv_layer updn'
!      write(*,*)mb,src,stype,rcv,updn
   if (n_layers_new.gt.nlay .OR. src.EQ.rcv) then
      write(0,*) n_layers_new,'>',nlay,' or source receiver at same level'
      call exit(1)
   endif
   nCom = 3 + 3*stype
   idx0 = 47
   if (stype.EQ.0) idx0 = 96

   flip = 1
   if ( rcv.GT.src ) then ! flip the model so rcv is above src
      flip = -1
      src = n_layers_new - src + 2
      rcv = n_layers_new - rcv + 2
   endif

   hs = 0.d0
!
   if(flip.lt.0)then
      call flip_model()
   endif
   call extra(hs)
   call extra2(hs)
!
   vs = new_vel_s(rcv) 
   call source(stype, xi(src), mu(src), flip)

! input sigma, number of points, dt, taper, and samples before first-arr
! sigma is the small imaginary part of the complex frequency, input in
! 1/time_length. Its role is to damp the trace (at rate of exp(-sigma*t))
! to reduce the wrape-arround.
!      write(0,'(a)') 'Input sigma NFFT dt lp nb smooth_factor wc1 wc2'
!      read(*,*) sigma,nfft,dt,taper,tb,smth,wc1,wc2
   if ( nfft.EQ.1 ) then
      dynamic = .FALSE.
      nfft2 = 1
      wc1 = 1
   else
      nfft2 = nfft/2
   endif
!      write(0,*)'check: ', wc1, wc2, smth, taper
   if (nfft2*smth .GT. nt) then
      write(0,'(a)') 'Number of points exceeds allowed'
      call exit(1)
   endif
   dw = twopi/(nfft*dt)
   sigma = sigma*dw/twopi
   wc = nfft2*(1.d0-taper)
   if (wc .LT. 1) wc=1
   taper = pi/(nfft2-wc+1)
   if (wc2.GT.wc) wc2=wc
   if (wc1.GT.wc2) wc1=wc2

! Input phase velocity window, dk, and kmax at w=0.
! pmin and pmax are in 1/vs.
! dk is in pi/max(x,hs). Since J(kx) oscillates with period 2pi/x at
! large distance, we require dk < 0.5 to guarentee at least 4 samples per
! period.
! kmax is in 1/hs. Because the kernels decay with k at rate of exp(-k*hs) at
! w=0, we require kmax > 10 to make sure we have have summmed enough.
!      write(0,'(a)') 'Input pmin pmax dk kc'
!      read(*,*) pmin,pmax,dk,kc
   pmin=0.d0
   pmax=1.11d0
   dk=0.2d0
   kc=30.d0

   kc = kc/hs
   pmin = pmin/vs
   pmax = pmax/vs

! input distance ranges
!      write(0,'(a)') 'Input number of distance ranges to compute'
!      read(*,*) nx
!      write(*,*) nx
   if (nx.gt.ndis) then
      write(0,'(a)') 'Number of ranges exceeds the max. allowed'
      call exit(1)
   endif
   xmax = hs
   do ix=1,nx
!	write(0, '(a)') 'Input x t0 output_name (2f10.3,1x,a)'
!        read(*,'(2f10.3,1x,a)')x(ix),t0(ix),fout(ix)
      if (xmax .LT. x(ix)) xmax=x(ix)
   enddo
   dk = dk*pi/xmax   ! satisfies bouchon criteria
   if(dk.gt.0.03)dk=0.03d0
   const = dk/twopi

!***************************************************************
!*************** do wavenumber integration for each frequency
   z = pmax*nfft2*dw/kc
   k = sqrt(z*z+1)
   total = nfft2*(kc/dk)*0.5*(k+log(k+z)/z)
!   write(0,'(a3,f9.5,a8,f9.2,a8,i9)')'dk',dk,'kmax',kc,'N',total
   tenPerc = 0.1*total
   count = 0.
   kc = kc*kc
   data1(:) = cmplx(0.,0.)
   

!   write(0,*)' start F-K computation, iw-range:',wc1,wc2,wc,nfft2
!$omp parallel &
!$omp& private (j, omega, w, att, k, n, u, i, ix, z, ixx, aj0, aj1, aj2, nf1, nf2), & 
!$omp& private (filter, phi, l, aux) &
!$omp& shared (wc1, nfft2, sigma, n_layers_new, pmin, dk, kc, pmax, src, rcv, mu, nx, x) &
!$omp& shared (aj0s, aj1s, aj2s, summ, flip, count, tenPerc, const, wc, wc2, taper, dynamic, t0, nCom)
!$omp do schedule(guided)
   do j=wc1, nfft2                                                  ! start frequency loop
      summ(:, :, j) = cmplx(0.0, 0.0)
      omega = (j-1)*dw
      w = cmplx(omega,-sigma,double)                             ! complex frequency
      w2 = w/twopi
      att = log(w2)/pi + cmplx(0.d0,0.5d0,double)                        ! A&R, p182
      do i = 1, n_layers_new
         ka(i) = w/(new_vel_p(i)*(1.d0+att/qqp(i)))
         kb(i) = w/(new_vel_s(i)*(1.d0+att/qqs(i)))
         ka(i) = ka(i)*ka(i)
         kb(i) = kb(i)*kb(i)
      enddo
      aux(:, :) = 0.d0
!
! replace summation in single precision to double precision, which reduces numerical error
!
      k = omega*pmin + 0.5d0*dk
      n = (sqrt(kc+(pmax*omega)**2)-k)/dk                           ! kmax
      do i=1,n                                                      ! start k-loop
         call kernel(k, u)
         do ix=1,nx
            z = k*x(ix)
            ixx = int(x(ix)/d_step + 0.001) + 1
!            call besselFn(z, aj0,aj1,aj2)
            aj0 = aj0s(ixx,i)
            aj1 = aj1s(ixx,i)
            aj2 = aj2s(ixx,i)
! n=0
            aux(1,ix) = aux(1,ix) + u(1,1)*aj0*flip
            aux(2,ix) = aux(2,ix) - u(2,1)*aj1
            aux(3,ix) = aux(3,ix) - u(3,1)*aj1
!            summ(1,ix,j) = summ(1,ix,j) + u(1,1)*aj0*flip
!            summ(2,ix,j) = summ(2,ix,j) - u(2,1)*aj1
!            summ(3,ix,j) = summ(3,ix,j) - u(3,1)*aj1
! n=1
!
!  l'hopital of bessel functions for small z
!
            if(z.gt.1.0e-10)then
               nf1 =    (u(2,2)+u(3,2))*aj1/z ! n=1
               nf2 = 2.*(u(2,3)+u(3,3))*aj2/z ! n=2
            else
               nf1 =    (u(2,2)+u(3,2))/2.   ! n=1
               nf2 = 2.*(u(2,3)+u(3,3))*z/8. ! n=2
            endif
            aux(4,ix) = aux(4,ix) + u(1,2)*aj1*flip
            aux(5,ix) = aux(5,ix) + u(2,2)*aj0 - nf1
            aux(6,ix) = aux(6,ix) + u(3,2)*aj0 - nf1
            aux(7,ix) = aux(7,ix) + u(1,3)*aj2*flip
            aux(8,ix) = aux(8,ix) + u(2,3)*aj1 - nf2
            aux(9,ix) = aux(9,ix) + u(3,3)*aj1 - nf2
!            summ(4,ix,j) = summ(4,ix,j) + u(1,2)*aj1*flip
!            summ(5,ix,j) = summ(5,ix,j) + u(2,2)*aj0 - nf1
!            summ(6,ix,j) = summ(6,ix,j) + u(3,2)*aj0 - nf1
!            summ(7,ix,j) = summ(7,ix,j) + u(1,3)*aj2*flip
!            summ(8,ix,j) = summ(8,ix,j) + u(2,3)*aj1 - nf2
!            summ(9,ix,j) = summ(9,ix,j) + u(3,3)*aj1 - nf2
         enddo
         k = k+dk
         count=count+1
         if ( mod(count,tenPerc) .EQ. 0) then
            write(0,'(i4,a6)') int(100.*count/total)+1, '% done'
         endif
      enddo                                                 ! end of k-loop
      summ(:, :nx, j) = aux(:, :nx)
      filter = const
      if ( dynamic .AND. j.GT.wc ) then 
         filter = 0.5d0*(1.d0+cos((j-wc)*taper))*filter
      endif
      if ( dynamic .AND. j.LT.wc2 ) then
         filter = 0.5d0*(1.d0+cos((wc2-j)*pi/(wc2-wc1)))*filter
      endif
      do ix=1,nx
         phi = omega*t0(ix)
         att = filter*cmplx(cos(phi),sin(phi),double)
         do l=1,nCom
            summ(l,ix,j) = summ(l,ix,j)*att
         enddo
      enddo
   enddo                            ! end of freqency loop
!$omp end do
!$omp end parallel
!      write(0,'(i9,a40)') count,' 100% done, writing files ... '
      
!***************************************************************
!*************** do inverse fourier transform
   dt = dt/smth
   nfft = smth*nfft
   nfft3 = nfft/2
   dfac = exp(sigma*dt)
   do ix=1,nx
      if ( nfft2.EQ.1 ) then
         do l=1,8
            green(1,l,ix) = real(summ(new_old(l),ix,1))
         enddo
      else
         if (disp) then
            do l=1,8
               do j=1,nfft2
                  data1(j) = summ(new_old(l),ix,j)
               enddo
               do j=nfft2+1,nfft3
                  data1(j) = cmplx(0.,0.)
               enddo
               call fftr(data1,nfft3,-dt)
               z2 = exp(sigma*t0(ix)) ! removing damping due sigma. Damping is w.r.t t=0
               do j=1,nfft3
                  tdata(2*j-1) = real(data1(j))*z2
                  z2 = z2*dfac
                  tdata(2*j) = aimag(data1(j))*z2
                  z2 = z2*dfac
               enddo
               green(1, l, ix)=0.0
               kahan_y(1, 1) = 0.0
               kahan_t(1, 1) = 0.0
               kahan_c(1, 1) = 0.0
               do j=2,nfft
                  kahan_y(1, 1) = (tdata(j) + tdata(j-1))*dt/2.0 - kahan_c(1, 1)
                  kahan_t(1, 1) = green(j-1,l,ix)+kahan_y(1, 1)
                  kahan_c(1, 1) = (kahan_t(1, 1)-green(j-1,l,ix))-kahan_y(1, 1)
                  green(j,l,ix) = kahan_t(1, 1)
               enddo
            enddo
         else
            do l=1,8
               do j=1,nfft2
                  data1(j) = summ(new_old(l),ix,j)
               enddo
               do j=nfft2+1,nfft3
                  data1(j) = cmplx(0.,0.)
               enddo
               call fftr(data1,nfft3,-dt)
               z2 = exp(sigma*t0(ix)) ! removing damping due sigma. Damping is w.r.t t=0
               do j=1,nfft3
                  tdata(2*j-1) = real(data1(j))*z2
                  z2 = z2*dfac
                  tdata(2*j) = aimag(data1(j))*z2
                  z2 = z2*dfac
               enddo
               green(1, l, ix)=0.0
               do j=1,nfft
                  green(j,l,ix)=tdata(j)
               enddo
            enddo
         endif
      endif
   enddo
   deallocate(summ)     
   end subroutine sub_bs_dc


end module fk_openmp
