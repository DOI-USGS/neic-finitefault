module vel_model_data


   use constants, only : nlay
   implicit none
   integer, private :: n_layers
   integer :: n_layers_new, s_layer
   real, private :: vel_p(nlay), vel_s(nlay), dens(nlay), thick(nlay), qp(nlay), qs(nlay)
   real :: new_vel_p(nlay), new_vel_s(nlay), new_dens(nlay), new_thick(nlay), qqp(nlay), qqs(nlay)
   real :: slowness(nlay)


contains


   subroutine read_vel_model(vel_model)
!
!  Args:
!  vel_model: file with velocity model
!
   implicit none
   character(len=100), intent(in) :: vel_model
   integer :: j
   open(1, file=vel_model, status='old')
   read(1, *) n_layers
   do j = 1, n_layers
      read (1, *)vel_p(j), vel_s(j), dens(j), thick(j), qp(j), qs(j)
   enddo
   end subroutine read_vel_model

  
   subroutine update_model(depth_source)
!
!  Args:
!  depth_source: depth of point source
!
   implicit none
   real, intent(in) :: depth_source
   real :: cum_dep, depths(nlay)
   integer :: j
   depths(1) = 0.0
   cum_dep = 0.0
   new_vel_p = vel_p
   new_vel_s = vel_s
   new_dens = dens
   new_thick = thick
   do j = 1, n_layers
      cum_dep = cum_dep + thick(j)
      depths(j + 1) = cum_dep
   enddo

   s_layer = 1
   do j = 1, n_layers
      if((depth_source .gt. depths(j)) .and. (depth_source .le. depths(j + 1))) then
         s_layer = j + 1
         exit
      end if
   enddo

   if (s_layer.gt.1) then
      new_thick(s_layer - 1) = depth_source - depths(s_layer - 1)
      new_thick(s_layer:n_layers + 1) = thick(s_layer - 1:n_layers)
      new_vel_p(s_layer:n_layers + 1) = vel_p(s_layer - 1:n_layers)
      new_vel_s(s_layer:n_layers + 1) = vel_s(s_layer - 1:n_layers)
      new_dens(s_layer:n_layers + 1) = dens(s_layer - 1:n_layers)
      new_thick(s_layer) = depths(s_layer) - depth_source
      if(new_thick(s_layer) .lt. 1.e-6)new_thick(s_layer) = 0.01
      n_layers_new = n_layers + 1
      s_layer = s_layer - 1
   endif

   slowness(:n_layers+1) = 1. / new_vel_p(:n_layers+1) ** 2

   end subroutine update_model


end module vel_model_data
