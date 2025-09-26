!     This program load the static subfault GFs and generate the matrix
!     for baseline correction.         Ji, Chen, 2003
!     Reviewed. The previous version have an error in loading GFs.
!                      
!                                   Ji, Chen, 2004
!  Fixed
!
!
!
module imagery_data


   use constants, only : max_subf, max_seg, max_subfaults, dpi
   implicit none
   integer, parameter, private :: max_points = 9000
   integer, private :: tracks, points
   integer :: ramp_length
   integer :: point_tracks(6)
   real :: weight_tracks(6)
   real*8 :: synm_whole(max_points), weight_sum
   real*8 :: ramp_gf(36, max_points)
   real :: lat(max_points), lon(max_points), max_los, syn_disp(max_points)
   real, private, allocatable :: green(:, :, :)
   real :: obse(max_points), look(3, max_points), weight(max_points)
   character(len=6) :: sta_name(max_points)
   character(len=10) :: ramp_types(6)
   integer :: subfaults, segments, nxs_sub(max_seg), nys_sub(max_seg)
  

contains


   subroutine get_imagery_gf()
   implicit none
   allocate(green(6, max_points, max_subfaults))
   end subroutine get_imagery_gf
  

   subroutine imagerydata_set_fault_parameters()
   use model_parameters, only : get_segments
   implicit none
   real :: dip(max_seg), strike(max_seg), delay_seg(max_seg)
   integer :: cum_subfaults(max_seg)
   call get_segments(nxs_sub, nys_sub, dip, strike, delay_seg, segments, subfaults, cum_subfaults)
   end subroutine imagerydata_set_fault_parameters
   

   subroutine get_imagery_data()
   implicit none
   integer i, no, k, points2, track
   integer cum_points(6)!point_tracks(6)
!   real weight_tracks(6)

   points2 = 0
   weight_sum = 0.0
   open(9, file='imagery_weights.txt', status='old')
   read(9,*) tracks
   do track = 1, tracks
      read(9,*) point_tracks(track), weight_tracks(track)
      weight_sum = weight_sum + point_tracks(track)*weight_tracks(track)
      points2 = points2 + point_tracks(track)
      cum_points(track) = points2
   enddo
   
!
   open(9, file='imagery_data.txt', status='old')
   read(9,*) points
   read(9,*)
   track = 1
   do i = 1,points
      read(9,*) no, sta_name(i), lat(i), lon(i), obse(i), (look(k, i), k = 1, 3)
      if (i .gt. cum_points(track)) track = track + 1
      weight(i) = weight_tracks(track)
   end do
   close(9) 
   max_los = maxval(abs(obse(:)))
   end subroutine get_imagery_data


   subroutine is_ramp(ramp_gf_file)
!
!  Args:
!  ramp_gf_file: True if there is a file with ramp GF data, False otherwise
!
   implicit none
   integer :: point, i, track
   logical :: ramp_gf_file
   character(len=10) :: ramp_type
   ramp_length = 0
   inquire( file = 'ramp_gf.txt', exist = ramp_gf_file )
   if (ramp_gf_file) then
      open(23, file='ramp_gf.txt', status='old')
      read(23,*) ramp_types(:tracks)
      do i = 1, tracks
         ramp_type = trim(ramp_types(i))
         select case (ramp_type)
         case ('static')
            ramp_length = ramp_length + 1
         case ('linear')
            ramp_length = ramp_length + 3
         case ('bilinear')
            ramp_length = ramp_length + 6
         case ('quadratic')
            ramp_length = ramp_length + 5
         end select
      enddo
      point = 0
      do track = 1, tracks
         ramp_type = trim(ramp_types(track))
         do i = 1, point_tracks(track)
            point = point + 1
            ramp_gf(:, point) = 0.d0
            if ((ramp_type .eq. 'None') .eqv. .False.) read(23,*)ramp_gf(:ramp_length, point)
         enddo
      enddo
      !do point = 1, points
      !   ramp_gf(:, point) = 0.d0
      !   read(23,*)ramp_gf(:ramp_length, point)
      !end do
      close(23)
   end if
   end subroutine is_ramp


   subroutine initial_ramp(ramp)
! 
!  Args: 
!  ramp: Value of imagery ramp
!
   implicit none
   real*8 :: ramp(36)
   integer :: i
   ramp(:) = 0.d0
   do i=1, ramp_length, 2
      ramp(i) = 20
   end do
   do i=2, ramp_length, 2
      ramp(i) = -20
   end do
   end subroutine initial_ramp


   subroutine initial_imagery(slip, rake, ramp)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  ramp: Value of imagery ramp, optional
!
   implicit none
   real*8, optional :: ramp(36)
   integer k, j, segment, point, subfault
   integer n_tt!, lines_asc, lines_desc
   real slip(:), rake(:)
   real :: cosal, sinal, angle
   real*8 :: disp, ramp2
   logical is_file
!
   open(33, file='Insar_static_subfault.txt', status='old')
   read(33,*) n_tt
   do point = 1, points
      read(33,*)
      subfault = 0
      do segment = 1, segments
         read(33,*)
         do j = 1, nys_sub(segment)*nxs_sub(segment)
            subfault = subfault + 1
            read(33,*)(green(k, point, subfault), k = 1, 6)
         enddo
      end do
   end do
   close(33)

   do point = 1, points
      syn_disp(point) = 0.d0
      do k = 1, 3
         j = 2 * k - 1
         disp = 0.d0
         do subfault = 1, subfaults
            angle = rake(subfault)*dpi 
            sinal = sin(angle)
            cosal = cos(angle)
            disp = disp + slip(subfault) &
     &  *(sinal*green(j, point, subfault)+cosal*green(j+1, point, subfault))
         end do
         disp = disp * look(k, point)
         syn_disp(point) = syn_disp(point) + disp
      end do
      if (present(ramp)) then
         ramp2 = 0.d0 
         do k = 1, ramp_length
            ramp2 = ramp2 + ramp_gf(k, point)*ramp(k)
         end do 
         syn_disp(point) = syn_disp(point) + ramp2
      endif
   end do
   open(10,file='imagery_synthetics.txt')
   write(10,*) points
   do point = 1, points
      write(10,*) point, sta_name(point), lat(point), lon(point), syn_disp(point)
   end do
   close(10)
   
   if (present(ramp)) then
      open(11,file='ramp_coefficients.txt')
      write(11,*) ramp(:)
      close(11)
   
      open(12,file='imagery_ramp.txt')
      write(12,*) points
      do point = 1, points
         ramp2 = 0.d0 
         do k = 1, ramp_length
            ramp2 = ramp2 + ramp_gf(k, point)*ramp(k)
         end do
         write(12,*) point, sta_name(point), lat(point), lon(point), ramp2
      end do
      close(12)
   end if
   end subroutine initial_imagery
  

!
! routine for loading static synthetic seismograms, given a rupture model
!
   subroutine imagery_synthetic(slip, rake, err, ramp)
!
!  Args:
!  slip: array with model slip values for all subfaults
!  rake: array with model rake values for all subfaults
!  err: imagery misfit
!  ramp: Value of imagery ramp, optional
!
   implicit none
   real*8, optional :: ramp(36)
   real slip(:), rake(:)!, err
   integer k, j, segment, point, subfault
   real err, dif, angle, sinal, cosal
   real*8 :: disp, err2, ramp2
      
   err2 = 0.d0
   synm_whole(:) = 0.d0
   do point = 1, points
      do k = 1, 3
         j = 2 * k - 1
         disp = 0.d0
         do subfault = 1, subfaults
            angle = rake(subfault)*dpi 
            sinal = sin(angle)
            cosal = cos(angle)
            disp = disp + slip(subfault) &
       &  *(sinal*green(j, point, subfault)+cosal*green(j+1, point, subfault))
         end do
         disp = disp * look(k, point)
         synm_whole(point) = synm_whole(point) + disp
      end do
      if (present(ramp)) then
         ramp2 = 0.d0 
         do k = 1, ramp_length
            ramp2 = ramp2 + ramp_gf(k, point)*ramp(k)
         end do 
         synm_whole(point) = synm_whole(point) + ramp2
      endif
      dif = synm_whole(point) - obse(point)
      err2 = err2 + weight(point) * dif * dif / max_los / max_los!100.0
   end do
   err2 = sqrt(err2/weight_sum)
   err = real(err2)
   end subroutine imagery_synthetic
            

!  
! subroutine for removing the static response of current subfault for all static stations
!  
   subroutine imagery_remove_subfault(slip, rake, subfault)
!
!  Args:
!  slip: slip value for this subfault
!  rake: rake value for this subfault
!  subfault: value of current subfault
!
   implicit none
   real, intent(in) :: slip, rake
   integer, intent(in) :: subfault
   integer k, j, point
   real disp, angle, sinal, cosal
!
   angle = rake*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   do point = 1, points
      disp = 0.0
      do k = 1, 3
         j = 2*k-1
         disp = disp + slip*look(k, point) &
     &   *(sinal*green(j, point, subfault)+cosal*green(j+1, point, subfault))
      end do
      synm_whole(point) = synm_whole(point)-disp
   end do

   end subroutine imagery_remove_subfault


!
! subroutine for testing the new response of the current subfault, for all stations
! we also give the misfit error of static data
!
   pure subroutine imagery_modify_subfault(slip, rake, subfault, err)
!
!  Args:
!  slip: slip value for this subfault
!  rake: rake value for this subfault
!  subfault: value of current subfault
!  err: imagery misfit
!
   implicit none
   real, intent(in) :: slip, rake
   real, intent(out) :: err
   integer, intent(in) :: subfault
   integer k, j, point
   real disp, angle, dif, sinal, cosal
   real*8 :: err2
!
   angle = rake*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   err2 = 0.d0
   do point = 1, points
      disp = 0.0
      do k = 1, 3
         j = 2*k-1
         disp = disp + slip*look(k, point) &
       & *(sinal*green(j, point, subfault)+cosal*green(j+1, point, subfault))
      end do
      dif = (synm_whole(point) + disp) - obse(point)
      err2 = err2 + weight(point) * dif * dif / max_los / max_los!100.0
   end do
   err2 = sqrt(err2/weight_sum)
   err = real(err2)

   end subroutine imagery_modify_subfault


!
! subroutine for asliping response of current subfault, for all stations.
! we also give the misfit error of static data
!
   subroutine imagery_add_subfault(slip, rake, subfault)
!
!  Args:
!  slip: slip value for this subfault
!  rake: rake value for this subfault
!  subfault: value of current subfault
!
   implicit none
   real, intent(in) :: slip, rake
   integer, intent(in) :: subfault
   integer k, j, point
   real disp, angle, dif, sinal, cosal

   angle = rake*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   do point = 1, points
      disp = 0.0
      do k = 1, 3
         j = 2*k-1
         disp = disp + look(k, point)*slip &
       & *(sinal*green(j, point, subfault)+cosal*green(j+1, point, subfault))
      end do
      synm_whole(point) = synm_whole(point)+disp
   end do

   end subroutine imagery_add_subfault


!  
! subroutine for removing the static response of current subfault for all static stations
!  
   subroutine imagery_remove_ramp(ramp)
!
!  Args:
!  ramp: Value of imagery ramp
!
   implicit none
   real*8, intent(in) :: ramp(36)
   integer k, j, point
   real ramp2
!
   do point = 1, points
      ramp2 = 0.0
      do k = 1, ramp_length
         ramp2 = ramp2 + ramp_gf(k, point)*ramp(k)
      end do
      synm_whole(point) = synm_whole(point)-ramp2
   end do

   end subroutine imagery_remove_ramp


!
! subroutine for testing the new response of the current subfault, for all stations
! we also give the misfit error of static data
!
   pure subroutine imagery_modify_ramp(ramp, err)
!
!  Args:
!  ramp: Value of imagery ramp
!  err: imagery misfit
!
   implicit none
   real*8, intent(in) :: ramp(36)
   real, intent(out) :: err
   integer k, j, point
   real*8 :: err2, ramp2, dif
!
   err2 = 0.d0
   do point = 1, points
      ramp2 = 0.0
      do k = 1, ramp_length
         ramp2 = ramp2 + ramp_gf(k, point)*ramp(k)
      end do
      dif = (synm_whole(point) + ramp2) - obse(point)
      err2 = err2 + weight(point) * dif * dif / max_los / max_los!100.0
   end do
   err2 = sqrt(err2/weight_sum)
   err = real(err2)

   end subroutine imagery_modify_ramp


!
! subroutine for asliping response of current subfault, for all stations.
! we also give the misfit error of static data
!
   subroutine imagery_add_ramp(ramp)
!
!  Args:
!  ramp: Value of imagery ramp
!
   implicit none
   real*8, intent(in) :: ramp(36)
   integer k, j, point
   real*8 :: err2, ramp2, dif
   
   do point = 1, points
      ramp2 = 0.d0
      do k = 1, ramp_length
         ramp2 = ramp2 + ramp_gf(k, point)*ramp(k)
      end do
      synm_whole(point) = synm_whole(point)+ramp2
   end do

   end subroutine imagery_add_ramp


   subroutine deallocate_imagery_gf()
   implicit none
   deallocate(green)
   end subroutine deallocate_imagery_gf


end module imagery_data
