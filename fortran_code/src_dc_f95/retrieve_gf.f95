!       Read Green's function from the Green Function bank 
!       Input dist_max dist_min d_step
!       Input dep_max dep_min dep_step
!       Input Lnpt dt
!       Input Green_name
!       output Green
!
module retrieve_gf


   implicit none
   integer :: lnpt, block_gg
   real :: dt, dep_max, dep_min, dep_step, dist_max, dist_min, d_step, t_cor


contains


   subroutine get_gf_data(gf_file, vel_model, gf_bank)
!
!  Args:
!  gf_file: file with properties of GF bank
!  vel_model: file with velocity model
!  gf_bank: location of GF bank
!
   implicit none
   character(len=100), intent(in) :: gf_file
   character(len=100), intent(out) :: vel_model, gf_bank
   open(1, file=gf_file, status='old')
   read(1, *)vel_model
   read(1, *)dep_max, dep_min, dep_step
   read(1, *)dist_max, dist_min, d_step
   read(1, *)lnpt, dt, block_gg, t_cor
   read(1, '(a)')gf_bank
   write(*,*)'GF bank name'
   write(*, *)gf_bank
   close(1)
   end subroutine get_gf_data


end module retrieve_gf
