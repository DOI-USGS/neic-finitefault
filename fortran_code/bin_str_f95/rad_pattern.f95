module rad_pattern


   use constants, only : dpi
   implicit none


contains

!
!     Eight Green´s function in the sequence SS DS SS DS DD SS DS DD
!
   subroutine rad_coef(dip, theta, az, ang_d, nn_comp, coef_v, coef_r)
!
!  Args:
!  dip: value of dip angle
!  theta: value of strike angle (?)
!  az: source-station azimuth
!  ang_d: angle to rotate radiation pattern
!  coef_v: vertical radiation coefficients
!  coef_r: radial radiation coefficients
!

   real*8 :: coef_v(2,3), coef_r(2,5), dip, theta, az, ang_d
   integer :: nn_comp
   real*8 :: source_az, source_az2, dip_rad, dip_rad2, rad_rot, tan_rot
   real*8 :: sind, cosd, sin2d, cos2d, sinaz, cosaz, sin2az, cos2az

   if(nn_comp.eq.16)then ! EW
      rad_rot = sin(ang_d)
      tan_rot = cos(ang_d)
   else ! NS
      rad_rot = cos(ang_d)
      tan_rot = -sin(ang_d)
   endif

   dip_rad = dip * dpi
   dip_rad2 = 2.d0 * dip_rad
   sind = sin(dip_rad)
   cosd = cos(dip_rad)
   sin2d = sin(dip_rad2)
   cos2d = cos(dip_rad2)

   source_az = (az - theta) * dpi
   source_az2 = 2.d0 * source_az
   sinaz = sin(source_az)
   cosaz = cos(source_az)
   sin2az = sin(source_az2)
   cos2az = cos(source_az2)
!
! vertical components
!
   coef_v(1, 1) = -0.5d0 * sin2d * cos2az ! vertical SS
   coef_v(1, 2) = -sinaz * cos2d ! vertical DS
   coef_v(1, 3) = 0.5d0 * sin2d ! vertical 45DS

   coef_v(2, 1) = -sind * sin2az ! vertical SS
   coef_v(2, 2) = cosaz * cosd ! vertical DS
   coef_v(2, 3) = 0.d0 ! vertical 45DS
!
! horizontal components
!
   coef_r(1, 1) = -0.5d0 * sin2d * sin2az * tan_rot ! tangential SS, rotated by rad_c
   coef_r(1, 2) = cosaz * cos2d * tan_rot ! tangential DS, rotated by rad_c
   coef_r(1, 3) = -0.5d0 * sin2d * cos2az * rad_rot ! radial SS, rotated by rad_c
   coef_r(1, 4) = -sinaz * cos2d * rad_rot ! radial DS, rotated by rad_c. 
   coef_r(1, 5) = 0.5d0 * sin2d * rad_rot ! radial 45DS rotated by rad_c

   coef_r(2, 1) = sind * cos2az * tan_rot ! tangential SS, rotated by rad_c
   coef_r(2, 2) = sinaz * cosd * tan_rot ! tangential DS, rotated by rad_c
   coef_r(2, 3) = -sind * sin2az * rad_rot ! radial SS, rotated by rad_c
   coef_r(2, 4) = cosaz * cosd * rad_rot ! radial DS, rotated by rad_c
   coef_r(2, 5) = 0.d0 ! radial 45, rotated by rad_c

   end subroutine rad_coef


end module rad_pattern
