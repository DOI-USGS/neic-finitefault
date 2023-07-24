module wave_travel
      

   implicit none 


contains


   subroutine trav(dis, nx, tmin)
!
!  Args:
!  dis: source-station distances to compute travel time
!  nx: amount of destinations for which to compute travel time
!  tmin: travel time at each source-station distance
!
!============================================================
! calculate travel time for horizontal layered model
! it outputs both times for first arrival and direct arrival
!
!============================================================
   use vel_model_data, only : s_layer, n_layers_new
   implicit none
   real, intent(out) :: tmin(*)
   real, intent(in) :: dis(*)
   integer, intent(in) :: nx
   real*8 :: t, t0, td, x
   complex*16 :: pd, p0, p
   integer :: i, lmax

   do i = 1, nx
      x = dis(i)
      lmax = s_layer
      call find2(x, lmax, td, pd)
      t0 = td
      p0 = pd
      do lmax = s_layer + 1, n_layers_new - 1
         call find2(x, lmax, t, p)
         if (t.lt.t0) then
            t0 = t
            p0 = p
         endif
      enddo
!c        write(*,'(f8.2,f8.2,f8.2,f8.4,f8.4)')x,t0,td,real(p0),real(pd)
      tmin(i) = t0
   enddo
   end subroutine trav


   subroutine find2(x, lmax, t0, p0)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     solve d(tau)/dp = 0 for p0, tau=p*x+eta*z
!     input:  x distance range
!     output: p0, t0
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   use vel_model_data, only : s_layer, slowness
   implicit none
   integer, intent(in) :: lmax
   integer :: i
   real*8, intent(in) :: x
   real*8, intent(out) :: t0
   complex*16, intent(out) :: p0 
   real*8 :: aminm, dtdp0, pm
   complex*16 :: p1, p2
!c
!c i am afraid i should make these of double precision
!c
   aminm = 1.d-6

   p1 = cmplx(0., 1.0d-20, kind(1d0))    ! make p slightly above the real axis

   pm = 9999999.d0
   do i = 1, lmax
      pm = amin1(pm, slowness(lmax))
   enddo
   p2 = sqrt(pm) + p1
   p0 = 0.5d0 * (p1 + p2)          ! the actual p0 lies betw. p1 and p2
   do while ( real(p2 - p1) .gt. aminm )
      dtdp0 = function2(x, p0, lmax)
      if( abs(dtdp0).lt.aminm )exit
      if( dtdp0.gt.0. )then
         p1 = p0
      else
         p2 = p0
      endif
      p0 = 0.5d0 * (p1 + p2)
   enddo
   pm  = amin1(pm, slowness(lmax+1))
   pm = sqrt(pm)
   if ((lmax.gt.s_layer) .and. (pm.lt.real(p0))) p0 = pm
   t0 = function1(p0, x, lmax)
   end subroutine find2


   function function1(p, x, lmax) result(taup)
!
!  Args:
!  p: 
!  x: source-station distance
!  lmax: layer of point source
!
   use vel_model_data, only : slowness, new_thick, s_layer
!c define function function1(p) = p x + eta h
   implicit none
   integer :: lmax, i
   real*8 :: x
   complex*16 :: p, taup
   complex :: pp
   taup = p * x
   pp = p * p
   do i = 1, s_layer
      taup = taup + sqrt(slowness(i) - pp) * new_thick(i)
   enddo
   do i = s_layer + 1, lmax
      taup = taup + sqrt(slowness(i) - pp) * new_thick(i) * 2.d0
   enddo
   end function function1


   function function2(x, p, lmax) result(dtdp)
!
!  Args:
!  p: 
!  x: distance
!  lmax: layer of point source
!
   use vel_model_data, only : slowness, new_thick, s_layer
!c define function2 = d(tau)/dp
   implicit none
   integer :: lmax, j
   real*8 :: x
   complex*16 :: p, dtdp
   complex :: pp
   pp = p * p
   dtdp = 0.d0
   do j = 1, s_layer
      dtdp = dtdp - new_thick(j) / sqrt(slowness(j) - pp)
   enddo
   do j = s_layer + 1, lmax
      dtdp = dtdp - 2.d0 * new_thick(j) / sqrt(slowness(j) - pp)
   enddo
   dtdp = x + p * dtdp
   end function function2


end module wave_travel
