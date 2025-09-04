!     This program load the static subfault GFs and generate the matrix
!     for baseline correction.         Ji, Chen, 2003
!     Reviewed. The previous version have an error in loading GFs.
!                      
!                                   Ji, Chen, 2004
!  Fixed
!
!
!
module static_data


   use constants, only : nnxy, max_seg, dpi
   use model_parameters, only : nxs_sub, nys_sub, delay_seg, dip_seg, stk_seg, nx_p, ny_p, n_seg
   implicit none
   integer, parameter, private :: n_stations = 500
   integer :: n_chan
   real*8 :: synm_whole(n_stations, 3), weight_sum
   real :: lat(n_stations), lon(n_stations)
   real :: green(n_stations, 6, nnxy, max_seg), syn_disp(n_stations, 3)
   real :: obse(n_stations, 3), weight(n_stations, 3)
   character(len=6) :: sta_name(n_stations)
  

contains


   subroutine initial_gnss(dd, aa)
   !!
   !!   Load static GF and compute synthetic displacements from a fault model.
   !!
   integer k, j, i_s, ir, no, nxy, i, i_seg
   integer iys, ixs, n_tt
   real dd(nnxy, max_seg), aa(nnxy, max_seg)
   real :: cosal, sinal, angle
   real*8 :: disp
   logical is_file
!
   inquire( file = 'Readlp.static', exist = is_file )
   if (is_file) then 
      weight_sum = 0.d0
      open(9,file='Readlp.static')
      read(9,*) n_chan
      read(9,*)
      do i = 1, n_chan
         read(9,*) no, sta_name(i), lat(i), lon(i), (obse(i, k), k = 1, 3), (weight(i, k), k = 1, 3)
         do k = 1, 3
            weight_sum = weight_sum + weight(i, k)
         end do
      end do
      close(9) 

      open(33,file='Green_static_subfault')
      read(33,*) n_tt
      do ir = 1, n_chan
         read(33,*)
         do i_seg = 1, n_seg
            read(33,*)
            nxy = 0
            do iys = 1, nys_sub(i_seg)
               do ixs = 1, nxs_sub(i_seg)
                  nxy = nxy + 1
                  read(33,*)(green(ir, k, nxy, i_seg), k = 1, 6)
               end do
            end do
         end do
      end do
      close(33)
      
      do k = 1, 3
         j = 2 * k - 1
         do ir = 1, n_chan
            disp = 0.d0
            do i_s = 1, n_seg
               nxy = 0
               do iys = 1, nys_sub(i_s)
                  do ixs = 1, nxs_sub(i_s)
                     nxy = nxy + 1
                     angle = aa(nxy, i_s)*dpi 
                     sinal = sin(angle)
                     cosal = cos(angle)
                     disp = disp + dd(nxy, i_s) &
     &  *(sinal*green(ir, j, nxy, i_s)+cosal*green(ir, j+1, nxy, i_s))
                  end do
               end do
            end do
            syn_disp(ir, k) = disp
         end do
      end do 
      open(10,file='synm.static')
      write(10,*) n_chan
      do ir = 1, n_chan
         write(10,*) ir, sta_name(ir), lat(ir), lon(ir),(syn_disp(ir, k), k=1, 3)
      end do
      close(10)
   end if
   end subroutine initial_gnss
   
   
   subroutine static_synthetic(dd, aa, nxys, err)
   !!
   !!   Load static synthetic displacements from a fault model.
   !!
   real dd(nnxy, max_seg), aa(nnxy, max_seg)!, err
   integer k, j, i_s, ir, nxys(max_seg), nxy
   real err, dif, angle, sinal, cosal
   real*8 :: disp, err2
      
   err2 = 0.d0
   synm_whole(:, :) = 0.d0
   do k = 1, 3
      j = 2 * k - 1
      do ir = 1, n_chan
         disp = 0.d0
         do i_s = 1, n_seg
            do nxy = 1, nxys(i_s)
               angle = aa(nxy, i_s)*dpi 
               sinal = sin(angle)
               cosal = cos(angle)
               disp = disp + dd(nxy, i_s) &
       &  *(sinal*green(ir, j, nxy, i_s)+cosal*green(ir, j+1, nxy, i_s))
            end do
         end do
         synm_whole(ir, k) = disp
         dif = synm_whole(ir, k) - obse(ir, k)
         err2 = err2 + weight(ir, k) * dif * dif / 100.0
      end do
   end do
   err2 = sqrt((err2/weight_sum))
   err = real(err2)
   end subroutine static_synthetic
            

   subroutine static_remove_subfault(dd, aa, n_s, n_sub)
   !!
   !!   Remove static response of current subfault, for all static stations.
   !!
   real, intent(in) :: dd(nnxy, max_seg), aa(nnxy, max_seg)
   integer, intent(in) :: n_s, n_sub
   integer k, j, ir
   real disp, angle, sinal, cosal
!
   angle = aa(n_sub, n_s)*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   do k = 1, 3
      j = 2*k-1
      do ir = 1, n_chan
         disp = dd(n_sub, n_s) &
     &   *(sinal*green(ir, j, n_sub, n_s)+cosal*green(ir, j+1, n_sub, n_s))
         synm_whole(ir, k) = synm_whole(ir, k)-disp
      end do
   end do

   end subroutine static_remove_subfault


   pure subroutine static_modify_subfault(slip, rake, n_s, n_sub, err)
   !!
   !!   We add the static response of the current subfault, for all static stations.
   !!   We compute misfit of static data
   !!
   real, intent(in) :: slip, rake
   real, intent(out) :: err
   integer, intent(in) :: n_s, n_sub
   integer k, j, ir
   real disp, angle, dif, sinal, cosal
   real*8 :: err2
!
   angle = rake*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   err2 = 0.d0
   do k = 1, 3
      j = 2*k-1
      do ir = 1, n_chan
         disp = slip &
       & *(sinal*green(ir, j, n_sub, n_s)+cosal*green(ir, j+1, n_sub, n_s))
         dif = (synm_whole(ir, k) + disp) - obse(ir, k)
         err2 = err2 + weight(ir, k) * dif * dif / 100.0
      end do
   end do
   err2 = sqrt(err2/weight_sum)
   err = real(err2)

   end subroutine static_modify_subfault


   subroutine static_add_subfault(dd, aa, n_s, n_sub, err)
   !!
   !!   Add static response of current subfault, for all stations.
   !!   Compute the misfit error of static data
   !!
   real, intent(in) :: dd(nnxy, max_seg), aa(nnxy, max_seg)
   real, intent(out) :: err
   integer, intent(in) :: n_s, n_sub
   integer k, j, ir
   real disp, angle, dif, sinal, cosal
   real*8 :: err2
!
   angle = aa(n_sub, n_s)*dpi
   sinal = sin(angle)
   cosal = cos(angle)

   err2 = 0.d0
   do k = 1, 3
      j = 2*k-1
      do ir = 1, n_chan
         disp = dd(n_sub, n_s) &
       & *(sinal*green(ir, j, n_sub, n_s)+cosal*green(ir, j+1, n_sub, n_s))
         synm_whole(ir, k) = synm_whole(ir, k)+disp
         dif = synm_whole(ir, k) - obse(ir, k)
         err2 = err2 + weight(ir, k) * dif * dif / 100.0
      end do
   end do
   err2 = sqrt(err2/weight_sum)
   err = real(err2)

   end subroutine static_add_subfault


end module static_data
