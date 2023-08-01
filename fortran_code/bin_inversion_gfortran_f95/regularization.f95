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


   subroutine slip_laplace(laplacian)
!
!  Args:
!  laplacian: laplacian slip regularization
!
!   Laplacian regularization of slip vector field
!
   implicit none
   real, intent(out) :: laplacian
   integer subfault, subfault2
   integer segment2, i, j
   real slip11, slip21, slip31, slip41!, error
   real slip12, slip22, slip32, slip42!, error
   real*8 :: laplace2, diff
   
   laplace2 = 0.d0
   do subfault = 1, subfaults
!                 write(*,*)"segment", segment, nx, ny, n_sub, jj
!       left
      slip11 = 0.0
      slip12 = 0.0
      slip21 = 0.0
      slip22 = 0.0
      slip31 = 0.0
      slip32 = 0.0
      slip41 = 0.0
      slip42 = 0.0
      segment2 = nleft(1, subfault)
      i = nleft(2, subfault)
      j = nleft(3, subfault)
      subfault2 = cum_subfaults(segment2)
      subfault2 = subfault2 + i+(j-1)*nxs_sub(segment2)
      if (segment2 .gt. 0) then
         slip11 = slip_field(1, subfault2)
         slip12 = slip_field(2, subfault2)
      end if
!        write(*,*) segment2, nxx, nyy, ll,"left"
!       right   
      segment2 = nright(1, subfault)
      i = nright(2, subfault)
      j = nright(3, subfault)
      subfault2 = cum_subfaults(segment2)
      subfault2 = subfault2 + i+(j-1)*nxs_sub(segment2)
      if (segment2 .gt. 0) then
         slip31 = slip_field(1, subfault2)
         slip32 = slip_field(2, subfault2)
      end if
!        write(*,*) segment2, nxx, nyy, ll,"right"
!       up    
      segment2 = nup(1, subfault)
      i = nup(2, subfault)
      j = nup(3, subfault)
      subfault2 = cum_subfaults(segment2)
      subfault2 = subfault2 + i+(j-1)*nxs_sub(segment2)
      if (segment2 .gt. 0) then
         slip21 = slip_field(1, subfault2)
         slip22 = slip_field(2, subfault2)
      end if
!                    write(*,*) segment2, nxx, nyy, ll,"up"
!       down
      segment2 = ndown(1, subfault)
      i = ndown(2, subfault)
      j = ndown(3, subfault)
      subfault2 = cum_subfaults(segment2)
      subfault2 = subfault2 + i+(j-1)*nxs_sub(segment2)
      if (segment2 .gt. 0) then
         slip41 = slip_field(1, subfault2)
         slip42 = slip_field(2, subfault2)
      end if
!                   write(*,*) segment2, nxx, nyy, ll,"down"
!
      diff = slip_field(1, subfault)-0.25*(slip11+slip21+slip31+slip41)
      diff = diff*diff
      laplace2 = laplace2+diff
      diff = slip_field(2, subfault)-0.25*(slip12+slip22+slip32+slip42)
      diff = diff*diff
      laplace2 = laplace2+diff
   end do
   laplace2 = sqrt(laplace2/subfaults)
   laplacian = real(laplace2)
   end subroutine slip_laplace


   pure subroutine time_laplace(rupt_time, laplacian)
!
!  Args:
!  rupt_time: rupture time at each subfault
!  laplacian: laplacian time regularization
!
!   Laplacian regularization of rupture initiation time
!
   implicit none
   real, intent(in) :: rupt_time(max_subfaults) 
   real, intent(out) :: laplacian
   integer segment2
   integer i, j
   real time1, time2, time3, time4!, error
   real*8 :: laplace2, diff
   integer subfault, subfault2

   laplace2 = 0.d0
   do subfault = 1, subfaults
      time1 = rupt_time(subfault)
      time2 = rupt_time(subfault)
      time3 = rupt_time(subfault)
      time4 = rupt_time(subfault)

!       left
      segment2 = nleft(1, subfault)
      i = nleft(2, subfault)
      j = nleft(3, subfault)
      subfault2 = cum_subfaults(segment2)
      subfault2 = subfault2 + i+(j-1)*nxs_sub(segment2)
      if (segment2 .gt. 0) then
         time1 = rupt_time(subfault2)
      end if
!       right   
      segment2 = nright(1, subfault)
      i = nright(2, subfault)
      j = nright(3, subfault)
      subfault2 = cum_subfaults(segment2)
      subfault2 = subfault2 + i+(j-1)*nxs_sub(segment2)
      if (segment2 .gt. 0) then
         time3 = rupt_time(subfault2)
      end if
!       up    
      segment2 = nup(1, subfault)
      i = nup(2, subfault)
      j = nup(3, subfault)
      subfault2 = cum_subfaults(segment2)
      subfault2 = subfault2 + i+(j-1)*nxs_sub(segment2)
      if (segment2 .gt. 0) then
         time2 = rupt_time(subfault2)
      end if
!       down
      segment2 = ndown(1, subfault)
      i = ndown(2, subfault)
      j = ndown(3, subfault)
      subfault2 = cum_subfaults(segment2)
      subfault2 = subfault2 + i+(j-1)*nxs_sub(segment2)
      if (segment2 .gt. 0) then
         time4 = rupt_time(subfault2)
      end if
      diff = rupt_time(subfault)-0.25*(time1+time2+time3+time4)
      diff = diff*diff
      laplace2 = laplace2+diff
   end do
   laplace2 = sqrt(laplace2/subfaults)
   laplacian = real(laplace2)
   end subroutine time_laplace


end module regularization
   
