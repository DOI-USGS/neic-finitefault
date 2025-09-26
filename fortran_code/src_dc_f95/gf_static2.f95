!
!	GNSS, we need use geometry coordanation directly.
!
program gf_static


   use constants, only : dpi, max_seg, nnxs, nnys, nnpx, nnpy, nt 
   use vel_model_data, only : read_vel_model, update_model, src, new_vel_s, new_dens
   use wave_travel, only : trav_fk
   use fk, only : sub_bs_dc 
   use rad_pattern, only : rad_coef
   use retrieve_gf, only : lnpt, dt
   use geodesics, only : distaz
   use bessel2, only : load_bessel
   implicit none
   character(len=10) sta_name
   character(len=9) vel_model
   real*8 :: ang_d(3000), dip, theta, azi(3000), bazi(3000)
   real lat_s(300),lon_s(300), area
   real dist(3000),lat_sta,lon_sta,lat_p,lon_p
   real dep_p,t0(3000)
   real az,baz,dis,ta0,dta
   real*8 :: coef_v(2,3),coef_r(2,5)
   real dxs,dys,niu
   integer j,dist_max,ir_max,ir,no,ix,iy
   integer iys,ixs,iyp,ixp,nxp,nyp,k,msou
   real green_p(nt,8,3000)
   real green_sub_dip(3), gf_dip
   real green_sub_stk(3), gf_stk
   real :: fau_mod(7, nnpx, nnpy, nnxs, nnys, max_seg)
   real :: green_dip(3, nnxs, nnys, max_seg, 300)
   real :: green_stk(3, nnxs, nnys, max_seg, 300)
   integer nxs0, nys0, n_seg,i_seg,ll
   integer nxs_sub(max_seg),nys_sub(max_seg)
   real :: dip_sub(max_seg),stk_sub(max_seg)
   real c_depth
   integer io_seg
   logical :: disp
   disp = .False.
!
! input position of gnss stations
!
!	pause
   open(12,file='Readlp.static')
   read(12,*)ir_max
   read(12,*)
   do ir=1,ir_max
      read(12,*)no,sta_name,lat_s(ir),lon_s(ir)
   enddo
   close(12)
!
!	Input the Latitude and Longtitude and depth of epicenter
!
   open(22,file='Fault.time')
   read(22,*)nxs0,nys0,c_depth
   read(22,*)n_seg,dxs,dys,nxp,nyp
   read(22,*)ta0,dta,msou
   do i_seg=1,n_seg
      write(*,*)ta0,dta 
      read(22,*)io_seg, dip_sub(i_seg),stk_sub(i_seg)
      read(22,*)nxs_sub(i_seg),nys_sub(i_seg)
      do iys=1,nys_sub(i_seg)
         do ixs=1,nxs_sub(i_seg)
            read(22,*)
         enddo
      enddo
   enddo
   close(22)
!
!	Input the position of point source and epicenter position
!
   open(22,file='Fault.pos')
   do i_seg=1,n_seg
      read(22,*)io_seg
      do iys=1,nys_sub(i_seg)
         do ixs=1,nxs_sub(i_seg)
            do iy=1,nyp
               do ix=1,nxp
                  read(22,*)(fau_mod(k,ix,iy,ixs,iys,i_seg),k=1,7)
               enddo
            enddo
         enddo
      enddo
   enddo
   close(22)
     
   lnpt = 0
   dt = 1.0
   area = dxs * dys
!
! input velocity model
!
!	open(20,file='debug')
   vel_model = 'vel_model'
   vel_model = trim(vel_model)
   call read_vel_model(vel_model)
   do ll=1,3000
      t0(ll)=-50
   enddo
   dist_max = 1001
!       
!       Compute static GF in all point sources
!
   do ir=1,ir_max
      lat_sta = lat_s(ir)
      lon_sta = lon_s(ir)
      do i_seg=1,n_seg
         dip=dip_sub(i_seg)
         theta=stk_sub(i_seg)
         write(*,*)'dip, strike', dip, theta
         do iys=1,nys_sub(i_seg)
            write(0,*)iys
            do iyp = 1, nyp
               dep_p = fau_mod(3,1,iyp,1,iys,i_seg)
               call update_model(dep_p,0.0)
               niu=new_vel_s(src)*new_vel_s(src)*new_dens(src)
!
!	   ++++++++++++++++++++++++++++++++++++++++++++++++
!
               ll = 0
               do ixs = 1, nxs_sub(i_seg)
                  do ixp = 1, nxp
                     ll = ll + 1
                     lat_p = fau_mod(1,ixp,iyp,ixs,iys,i_seg)
                     lon_p = fau_mod(2,ixp,iyp,ixs,iys,i_seg)
                     call distaz(lat_sta,lon_sta,lat_p,lon_p,dis,az,baz)
                     dist(ll)=dis
                     azi(ll)=az
                     bazi(ll)=baz
                     dist(ll)=max(0.25, dist(ll))
                     ang_d(ll)=(baz - 180.0)*dpi
                  enddo
               enddo
               call sub_bs_dc(ll,dist,t0,green_p,disp)
               ll = 0
               green_dip(:,:,iys,i_seg,ir) = 0.0
               green_stk(:,:,iys,i_seg,ir) = 0.0
               do ixs = 1, nxs_sub(i_seg)
                  do ixp = 1, nxp
                     ll = ll + 1
                     gf_dip = 0.0
                     gf_stk = 0.0
!
! Vertical component
!             
                     call rad_coef(dip,theta,azi(ll),ang_d(ll),14,coef_v,coef_r)
                     do j=1,3
                        gf_dip = gf_dip + coef_v(1,j)*green_p(1,j+5,ll)*niu*area
                        gf_stk = gf_stk + coef_v(2,j)*green_p(1,j+5,ll)*niu*area
                     enddo
                     green_dip(1,ixs,iys,i_seg,ir)=gf_dip + green_dip(1,ixs,iys,i_seg,ir)
                     green_stk(1,ixs,iys,i_seg,ir)=gf_stk + green_stk(1,ixs,iys,i_seg,ir)
                     gf_dip = 0.0
                     gf_stk = 0.0
!
! horizontal components
!
                     call rad_coef(dip,theta,azi(ll),ang_d(ll),15,coef_v,coef_r) !N
                     do j=1,5
                        gf_dip = gf_dip + coef_r(1,j)*green_p(1,j,ll)*niu*area       
                        gf_stk = gf_stk + coef_r(2,j)*green_p(1,j,ll)*niu*area
                     enddo
                     green_dip(2,ixs,iys,i_seg,ir)=gf_dip + green_dip(2,ixs,iys,i_seg,ir)
                     green_stk(2,ixs,iys,i_seg,ir)=gf_stk + green_stk(2,ixs,iys,i_seg,ir)
                     gf_dip = 0.0
                     gf_stk = 0.0
                     call rad_coef(dip,theta,azi(ll),ang_d(ll),16,coef_v,coef_r) !E
                     do j=1,5
                        gf_dip = gf_dip + coef_r(1,j)*green_p(1,j,ll)*niu*area  
                        gf_stk = gf_stk + coef_r(2,j)*green_p(1,j,ll)*niu*area
                     enddo
                     green_dip(3,ixs,iys,i_seg,ir)=gf_dip + green_dip(3,ixs,iys,i_seg,ir)
                     green_stk(3,ixs,iys,i_seg,ir)=gf_stk + green_stk(3,ixs,iys,i_seg,ir)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo

!
!       Save Green functions for each subfault
!
   open(13,file='Green_static_subfault')
   write(13,*)ir_max
   do ir=1,ir_max
      write(13,*)ir
      do i_seg=1,n_seg
         write(13,*)i_seg
         do iys=1,nys_sub(i_seg)
            do ixs=1,nxs_sub(i_seg)
               do j = 1, 3
                  green_sub_dip(j) = green_dip(j,ixs,iys,i_seg,ir) / float(nxp*nyp)
                  green_sub_stk(j) = green_stk(j,ixs,iys,i_seg,ir) / float(nxp*nyp)
               enddo
               write(13,*)green_sub_dip(1),green_sub_stk(1),green_sub_dip(2), &
                       &  green_sub_stk(2),green_sub_dip(3),green_sub_stk(3)
            enddo
         enddo
      enddo
   enddo
      
   close(13)


end program gf_static
