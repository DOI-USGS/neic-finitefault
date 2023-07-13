module bessel2


   use constants, only : pi, twopi, ndis
   use retrieve_gf, only : dist_max, lnpt, dt
   implicit none
   integer, parameter :: dmax=2000, fmax=500000
   real, allocatable :: aj0s(:, :), aj1s(:, :), aj2s(:, :)


contains


   subroutine load_bessel(x,nd_max)
!
!  Args:
!  x: range of source-station distances
!  nd_max: amount of source-station distances
!
   implicit none
   real :: aj0, aj1, aj2, x(ndis)
   real :: k, dk, z, pmax, omega, dw
   integer :: i, n, ix, kc, nfft2, nfft, nd_max
   allocate(aj0s(dmax, fmax))
   allocate(aj1s(dmax, fmax))
   allocate(aj2s(dmax, fmax))
   write(0,*)'Load bessel functions to memory...'
   dk = 0.2
   dk = dk*pi/dist_max
   if(dk.gt.0.03)dk=0.03
   k=0.5*dk 
   kc = 30 / 0.5
   nfft = 2 ** lnpt
   nfft2 = nfft
   dw = twopi/(nfft*dt)
   omega=(nfft2-1)*dw
   pmax = 1.11
   n=(sqrt(kc*kc+(pmax*omega)**2)-k)/dk
   write(*,*)"maximum n is ",n
   if(n.gt.fmax)then
      write(*,*)n,fmax
      write(*,*)"You might want to increase the size of nkk_max"
      n=fmax
   endif
   do i=1,n
      do ix=1,nd_max
         z=x(ix)*k
         call besselFn(z,aj0,aj1,aj2)
         aj0s(ix,i)=aj0
         aj1s(ix,i)=aj1
         aj2s(ix,i)=aj2
      enddo
      k=k+dk
   enddo

   end subroutine load_bessel

   subroutine deallocate_bessel()
   deallocate(aj0s)
   deallocate(aj1s)
   deallocate(aj2s)
   end subroutine deallocate_bessel 


end module bessel2
