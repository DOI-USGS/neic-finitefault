module regularization
   

   use constants, only : max_seg, max_subf, max_subfaults, dpi
   implicit none
   real :: slip_field(2, max_subfaults)
   real :: rake_min
   integer :: subfaults, cum_subfaults(max_seg)
   integer :: nxs_sub(max_seg), nys_sub(max_seg)
   integer :: nleft(3, max_subfaults), nright(3, max_subfaults), & 
   & nup(3, max_subfaults), ndown(3, max_subfaults)


contains


   subroutine regularization_set_fault_parameters()
   use model_parameters, only : get_segments, get_borders
   implicit none
   real :: dip(max_seg), strike(max_seg), delay_seg(max_seg)
   integer :: segments
   call get_segments(nxs_sub, nys_sub, dip, strike, delay_seg, segments, subfaults, cum_subfaults)
   call get_borders(rake_min, nleft, nright, nup, ndown)   
   end subroutine regularization_set_fault_parameters
   
        
   subroutine define_slip_field(slip, rake)
!
!  Args:
!  slip: slip at each subfault
!  rake: rake at each subfault
!
!   Laplacian regularization of slip vector field
!
   implicit none
   real, intent(in) :: slip(max_subfaults), rake(max_subfaults)
   real angle
   integer subfault
!  
!  we define the slip vector field
!  
   do subfault = 1, subfaults
      angle = (rake(subfault)-rake_min)*dpi
      slip_field(1, subfault) = slip(subfault)*cos(angle)
      slip_field(2, subfault) = slip(subfault)*sin(angle)
   end do
   end subroutine define_slip_field

   
   subroutine modify_slip_field(subfault, slip_subf, rake_subf)
!
!  Args:
!  subfault: current subfault
!  slip_subf: slip at current subfault
!  rake_subf: rake at current subfault
!
!   Laplacian regularization of slip vector field
!
   implicit none
   integer, intent(in) :: subfault
   real, intent(in) :: slip_subf, rake_subf
   real angle
   angle = (rake_subf-rake_min)*dpi
   slip_field(1, subfault) = slip_subf*cos(angle)
   slip_field(2, subfault) = slip_subf*sin(angle)
   end subroutine modify_slip_field


   subroutine slip_laplace(err)
!
!  Args:
!  err: laplacian slip regularization
!
!   Laplacian regularization of slip vector field
!
   implicit none
   real, intent(out) :: err
   integer subfault, subfault2
   integer n_is, nxx, nyy
   real d11, d21, d31, d41!, error
   real d12, d22, d32, d42!, error
   real*8 :: err2, error, err3, err4
   
   err2 = 0.d0
   do subfault = 1, subfaults
!                 write(*,*)"segment", segment, nx, ny, n_sub, jj
!       left
      d11 = 0.0
      d12 = 0.0
      d21 = 0.0
      d22 = 0.0
      d31 = 0.0
      d32 = 0.0
      d41 = 0.0
      d42 = 0.0
      n_is = nleft(1, subfault)
      nxx = nleft(2, subfault)
      nyy = nleft(3, subfault)
      subfault2 = cum_subfaults(n_is)
      subfault2 = subfault2 + nxx+(nyy-1)*nxs_sub(n_is)
      if (n_is .gt. 0) then
         d11 = slip_field(1, subfault2)
         d12 = slip_field(2, subfault2)
      end if
!        write(*,*) n_is, nxx, nyy, ll,"left"
!       right   
      n_is = nright(1, subfault)
      nxx = nright(2, subfault)
      nyy = nright(3, subfault)
      subfault2 = cum_subfaults(n_is)
      subfault2 = subfault2 + nxx+(nyy-1)*nxs_sub(n_is)
      if (n_is .gt. 0) then
         d31 = slip_field(1, subfault2)
         d32 = slip_field(2, subfault2)
      end if
!        write(*,*) n_is, nxx, nyy, ll,"right"
!       up    
      n_is = nup(1, subfault)
      nxx = nup(2, subfault)
      nyy = nup(3, subfault)
      subfault2 = cum_subfaults(n_is)
      subfault2 = subfault2 + nxx+(nyy-1)*nxs_sub(n_is)
      if (n_is .gt. 0) then
         d21 = slip_field(1, subfault2)
         d22 = slip_field(2, subfault2)
      end if
!                    write(*,*) n_is, nxx, nyy, ll,"up"
!       down
      n_is = ndown(1, subfault)
      nxx = ndown(2, subfault)
      nyy = ndown(3, subfault)
      subfault2 = cum_subfaults(n_is)
      subfault2 = subfault2 + nxx+(nyy-1)*nxs_sub(n_is)
      if (n_is .gt. 0) then
         d41 = slip_field(1, subfault2)
         d42 = slip_field(2, subfault2)
      end if
!                   write(*,*) n_is, nxx, nyy, ll,"down"
!
      error = slip_field(1, subfault)-0.25*(d11+d21+d31+d41)
      error = error*error
      err2 = err2+error
      error = slip_field(2, subfault)-0.25*(d12+d22+d32+d42)
      error = error*error
      err2 = err2+error
   end do
   err2 = sqrt(err2/subfaults)
   err = real(err2)
   end subroutine slip_laplace


   pure subroutine time_laplace(tt, err)
!
!  Args:
!  tt: rupture time at each subfault
!  err: laplacian time regularization
!
!   Laplacian regularization of rupture initiation time
!
   implicit none
   real, intent(in) :: tt(max_subfaults) 
   real, intent(out) :: err
   integer n_is
   integer nxx, nyy
   real d1, d2, d3, d4!, error
   real*8 :: err2, error
   integer subfault, subfault2, j

   err2 = 0.d0
   do subfault = 1, subfaults
      d1 = tt(subfault)
      d2 = tt(subfault)
      d3 = tt(subfault)
      d4 = tt(subfault)

!       left
      n_is = nleft(1, subfault)
      nxx = nleft(2, subfault)
      nyy = nleft(3, subfault)
      subfault2 = cum_subfaults(n_is)
      subfault2 = subfault2 + nxx+(nyy-1)*nxs_sub(n_is)
      if (n_is .gt. 0) then
         d1 = tt(subfault2)
      end if
!       right   
      n_is = nright(1, subfault)
      nxx = nright(2, subfault)
      nyy = nright(3, subfault)
      subfault2 = cum_subfaults(n_is)
      subfault2 = subfault2 + nxx+(nyy-1)*nxs_sub(n_is)
      if (n_is .gt. 0) then
         d3 = tt(subfault2)
      end if
!       up    
      n_is = nup(1, subfault)
      nxx = nup(2, subfault)
      nyy = nup(3, subfault)
      subfault2 = cum_subfaults(n_is)
      subfault2 = subfault2 + nxx+(nyy-1)*nxs_sub(n_is)
      if (n_is .gt. 0) then
         d2 = tt(subfault2)
      end if
!       down
      n_is = ndown(1, subfault)
      nxx = ndown(2, subfault)
      nyy = ndown(3, subfault)
      subfault2 = cum_subfaults(n_is)
      subfault2 = subfault2 + nxx+(nyy-1)*nxs_sub(n_is)
      if (n_is .gt. 0) then
         d4 = tt(subfault2)
      end if
      error = tt(subfault)-0.25*(d1+d2+d3+d4)
      error = error*error
      err2 = err2+error
   end do
   err2 = sqrt(err2/subfaults)
   err = real(err2)
   end subroutine time_laplace


end module regularization
   
