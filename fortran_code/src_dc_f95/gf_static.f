c
c	GNSS, we need use geometry coordanation directly.
c
	implicit none
        include 'omp_lib.h'
	include 'fk.h'
	character*10 sta_name
	real dpi
	parameter(dpi=3.1415926535898/180.d0)
	integer max_seg, nnpx, nnpy, nnxs, nnys, nlay_r, nlay_s
	parameter(max_seg=5, nnpx=15, nnpy=15, nnxs=30, nnys=20)
	real ang_d(300), lat_s(300),lon_s(300), area
	real dist(300),lat_sta,lon_sta,lat_p,lon_p
	real azi(300),bazi(300),dip,rake,theta
	real dep_p,dt,t0(300)
	real az,baz,dis,ta0,dta
	real coef_v(2,3),coef_r(2,5),x,rad_c
	real dxs,dys,v_min
	real t1,w,niu
	integer i,j,block_gg,nsta, dist_max
	integer ir_max, n_chan,ir,nx,ky,kx,no
	integer iys,ixs,iyp,ixp,nxp,nyp
	integer npt,lnpt,npxy,k
	integer msou
	integer io_v_d
	real green_p(nt,8,ndis)
	real green_s(nt,8),v_ref
c	real green_dip(3,300,nnxs,nnys,max_seg), gf_dip
c	real green_stk(3,300,nnxs,nnys,max_seg), gf_stk
	real green_sub_dip(3), gf_dip
	real green_sub_stk(3), gf_stk
        real :: fau_mod(7, nnpx, nnpy, nnxs, nnys, max_seg)
        real, allocatable :: fau_mod2(:, :, :, :)
        real :: green_dip(3, 300, nnxs, nnys, max_seg)
        real, allocatable :: green_dip2(:, :, :)
        real :: green_stk(3, 300, nnxs, nnys, max_seg)
        real, allocatable :: green_stk2(:, :, :)
c	dimension fau_mod(7,nnpx,nnpy,nnxs,nnys,max_seg)
	integer nxs0, nys0, n_seg,i_seg
	integer ix,iy,nxs_sub,nys_sub

	real dip_sub,stk_sub,delay_sub
	real c_depth,t_cor
	dimension dip_sub(Max_seg),stk_sub(Max_seg)
	dimension nxs_sub(Max_seg),nys_sub(Max_seg),delay_sub(Max_seg)	
	integer jo,jjo,nb,io_seg
	real c(0:nlay),s(0:nlay),den(0:nlay),th(0:nlay),qa(0:nlay),qb(0:nlay)
	real vel_p(0:nlay),vel_s(0:nlay),dens(0:nlay),thick(0:nlay),
     *  qqa(0:nlay),qqb(0:nlay)
	real dist_min, max_dist, dep_min, dep_max, z_min, z_max
	common/model_fk/jo,vel_p,vel_s,dens,thick, qqa, qqb

	real d_step, dep_step
	common/gr_fun_prop/npt,dt,t_cor        
	common/bounds/dist_min,dist_max,dep_min,dep_max,d_step,dep_step
        allocate(fau_mod2(7, nnpx, nnpy, nnxs))
        allocate(green_dip2(3, 300, nnxs))
        allocate(green_stk2(3, 300, nnxs))
c
c input position of gnss stations
c
c	pause
	open(12,file='Readlp.static')
	read(12,*)ir_max
	read(12,*)
	do ir=1,ir_max
	   read(12,*)no,sta_name,lat_s(ir),lon_s(ir)
	enddo
	close(12)
c
c	Input the Latitude and Longtitude and depth of epicenter
c
	open(22,file='Fault.time')
	read(22,*)NXS0,NYS0,c_depth
	read(22,*)n_seg,dxs,dys,nxp,nyp
	read(22,*)TA0,DTA,MSOU,v_ref,io_v_d
	do i_seg=1,n_seg
	   read(22,*)io_seg, dip_sub(i_seg),stk_sub(i_seg)
	   read(22,*)nxs_sub(i_seg),nys_sub(i_seg),delay_sub(i_seg)
	   do iys=1,nys_sub(i_seg)
	      do ixs=1,nxs_sub(i_seg)
	         read(22,*)
	      enddo
	enddo
	enddo
	close(22)
c
c	Input the position of point source and epicenter position
c
	open(22,file='Fault.pos')
	do i_seg=1,n_seg
	   read(22,*)io_seg
	   do 55 iys=1,nys_sub(i_seg)
	      do 56 ixs=1,nxs_sub(i_seg)
	         do 57 iy=1,nyp
	            do 58 ix=1,nxp
	               read(22,*)(fau_mod(k,ix,iy,ixs,iys,i_seg),k=1,7)
 58                 continue
 57              continue
 56           continue
 55        continue
	enddo
	close(22)
    
	lnpt = 0
	dt = 1.0
	area = dxs * dys
c
c input velocity model
c
c	open(20,file='debug')
	open(11,file='vel_model',status='old')
	read(11,*) jo
	write(*,*)'jo: ', jo
	do j=1,jo
	   read(11,*) vel_p(j),vel_s(j),dens(j),thick(j),qqa(j),qqb(j)
	   write(*,*) vel_p(j),vel_s(j),dens(j),thick(j),qqa(j),qqb(j)
	enddo
	if(thick(jo).gt.1.0e-10)then
	   jo=jo+1
	   thick(jo)=0.0
	   vel_p(jo)=vel_p(jo-1)
	   vel_s(jo)=vel_s(jo-1)
	   dens(jo)=dens(jo-1)
	   qqa(jo)=qqa(jo-1)
	   qqb(jo)=qqb(jo-1)
	endif
	close(11)
        do ir=1,300
           t0(ir)=-50
        enddo
        dist_max = 1001
c       
c       Compute static GF in all point sources
c
	do i_seg=1,n_seg
	   dip=dip_sub(i_seg)
	   theta=stk_sub(i_seg)
	   write(*,*)'dip, strike', dip, theta
	   do iys=1,nys_sub(i_seg)
              write(0,*)iys
              fau_mod2(:, :, :, :) = fau_mod(:, :, :, :, iys, i_seg)
c$omp parallel
c$omp& shared(fau_mod2, dip, theta, i_seg, nxs_sub, nyp, nxp, lat_s,
c$omp+ lon_s, ir_max, lnpt, dt, t0, dist_max, area)
c$omp& private(ixs, j, green_dip2, green_stk2, iyp, dep_p, jjo, nlay_s,
c$omp+ nlay_r, c, s, den, th, qa, qb, niu, ixp, lat_p, lon_p, ir,
c$omp+ lat_sta, lon_sta, dis, az, baz, dist, azi, bazi, ang_d, green_p,
c$omp+ gf_dip, gf_stk, coef_v, coef_r)
c$omp+ lon_s, ir_max, lnpt, dt, t0, dist_max, area)
c$omp do
	      do ixs=1,nxs_sub(i_seg)
c
c	   ++++++++++++++++++++++++++++++++++++++++++++++++
c
		 do j=1,3
		    green_dip2(j,ir,ixs)=0.0
		    green_stk2(j,ir,ixs)=0.0
		 enddo
		 do iyp = 1, nyp
		    dep_p = fau_mod2(3,1,iyp,ixs)
		    call cmodel_fk(dep_p,0.0,jjo,nlay_s,nlay_r,c,s,den,th,qa,qb)
		    niu=s(nlay_s)*s(nlay_s)*den(nlay_s)
		    do ixp = 1, nxp
		       lat_p=fau_mod2(1,ixp,iyp,ixs)
		       lon_p=fau_mod2(2,ixp,iyp,ixs)
		       do ir=1,ir_max
			  lat_sta = lat_s(ir)
			  lon_sta = lon_s(ir)
			  call distaz(lat_sta,lon_sta,lat_p,lon_p,dis,az,baz)
			  dist(ir)=dis
			  azi(ir)=az
			  bazi(ir)=baz
			  dist(ir)=max(0.25, dist(ir))
			  ang_d(ir)=(baz - 180.0)*dpi
		       enddo
		       call sub_bs_dc(ir_max,dist,lnpt,dt,t0,green_p,dist_max,dep_p)
		       do ir=1,ir_max
			  gf_dip = 0.0
			  gf_stk = 0.0
c
c Vertical component
c             
			  call rad_coef(dip,theta,azi(ir),ang_d(ir),14,coef_v,coef_r)
			  do j=1,3
			     gf_dip = gf_dip + coef_v(1,j)*green_p(1,j+5,ir)*niu*area
			     gf_stk = gf_stk + coef_v(2,j)*green_p(1,j+5,ir)*niu*area
			  enddo
			  green_dip2(1,ir,ixs)=gf_dip + green_dip2(1,ir,ixs)
			  green_stk2(1,ir,ixs)=gf_stk + green_stk2(1,ir,ixs)
			  gf_dip = 0.0
			  gf_stk = 0.0
c
c horizontal components
c
			  call rad_coef(dip,theta,azi(ir),ang_d(ir),15,coef_v,coef_r) !N
			  do j=1,5
			     gf_dip = gf_dip + coef_r(1,j)*green_p(1,j,ir)*niu*area       
			     gf_stk = gf_stk + coef_r(2,j)*green_p(1,j,ir)*niu*area
			  enddo
			  green_dip2(2,ir,ixs)=gf_dip + green_dip2(2,ir,ixs)
			  green_stk2(2,ir,ixs)=gf_stk + green_stk2(2,ir,ixs)
			  gf_dip = 0.0
			  gf_stk = 0.0
			  call rad_coef(dip,theta,azi(ir),ang_d(ir),16,coef_v,coef_r) !E
			  do j=1,5
			     gf_dip = gf_dip + coef_r(1,j)*green_p(1,j,ir)*niu*area  
			     gf_stk = gf_stk + coef_r(2,j)*green_p(1,j,ir)*niu*area
			  enddo
			  green_dip2(3,ir,ixs)=gf_dip + green_dip2(3,ir,ixs)
			  green_stk2(3,ir,ixs)=gf_stk + green_stk2(3,ir,ixs)
		       enddo
		    enddo
		 enddo
	      enddo
c$omp end do
c$omp end parallel
              green_dip(:, :, :, iys, i_seg) = green_dip2(:, :, :)
              green_stk(:, :, :, iys, i_seg) = green_stk2(:, :, :)
 	   enddo
        enddo

        deallocate(fau_mod2)
        deallocate(green_dip2)
        deallocate(green_stk2)
c
c       Save Green functions for each subfault
c
	open(13,file='Green_static_subfault')
	write(13,*)ir_max
	do ir=1,ir_max
	   write(13,*)ir
	   do i_seg=1,n_seg
	      write(13,*)i_seg
	      do iys=1,nys_sub(i_seg)
		 do ixs=1,nxs_sub(i_seg)
		    do j = 1, 3
		       green_sub_dip(j) = green_dip(j,ir,ixs,iys,i_seg) /
     *   float(nxp*nyp)
		       green_sub_stk(j) = green_stk(j,ir,ixs,iys,i_seg) /
     *   float(nxp*nyp)
		    enddo
		 write(13,*)
     *   green_sub_dip(1),green_sub_stk(1),green_sub_dip(2),
     *   green_sub_stk(2),green_sub_dip(3),green_sub_stk(3)
		 enddo
	      enddo
	   enddo
	enddo
      
	close(13)
	stop
	end


c
c     Eight Green´s function in the sequence SS DS SS DS DD SS DS DD
c
	subroutine rad_coef(dip,theta,az,ang_d,nn_comp,coef_v,coef_r)
c
	parameter(dpi=3.1415926535898/180.d0)
        real coef_v,coef_r,d,dip,theta,az,ang_d
        real aaz,cosaaz,sinaaz,cos2aaz,sin2aaz,cos2d,sin2d
        real cosd,sind,d2,rad_rot,tan_rot
        integer nn_comp
        dimension coef_v(2,3),coef_r(2,5)
c
c       Radiation pattern for strong motion GF
c
	 if(nn_comp.eq.16)then ! EW
	      rad_rot=sin(ang_d)
	      tan_rot=cos(ang_d)
	 elseif(nn_comp.le.15)then ! NS
	      rad_rot=cos(ang_d)
	      tan_rot=-sin(ang_d)
	 endif

         D=Dip*dpi
         COSD=COS(D)
         SIND=SIN(D)

         D2=2.0*D
         COS2D=COS(D2)
         SIN2D=SIN(D2)

         AAZ=(AZ-THETA)*dpi
         COSAAZ=COS(AAZ)
         SINAAZ=SIN(AAZ)

         AAZ=2.0*AAZ
         COS2AAZ=COS(AAZ)
         SIN2AAZ=SIN(AAZ)
c
c vertical components
c
         coef_v(1,1)=-0.5*sin2D*cos2AAZ ! vertical SS
         coef_v(1,2)=-sinAAZ*cos2D ! vertical DS
         coef_v(1,3)=0.5*sin2D ! vertical 45DS

         coef_v(2,1)=-sinD*sin2AAZ ! vertical SS
         coef_v(2,2)=cosAAZ*cosD ! vertical DS
         coef_v(2,3)=0 ! vertical 45DS
c
c horizontal components
c
         coef_r(1,1)=-0.5*sin2D*sin2AAZ*tan_rot ! tangential SS, rotated by rad_c
         coef_r(1,2)=cosAAZ*cos2D*tan_rot ! tangential DS, rotated by rad_c
         coef_r(1,3)=-0.5*sin2D*cos2AAZ*rad_rot ! radial SS, rotated by rad_c
         coef_r(1,4)=-sinAAZ*cos2D*rad_rot ! radial DS, rotated by rad_c. 
         coef_r(1,5)=0.5*sin2D*rad_rot ! radial 45DS rotated by rad_c

         coef_r(2,1)=sinD*cos2AAZ*tan_rot ! tangential SS, rotated by rad_c
         coef_r(2,2)=sinAAZ*cosD*tan_rot ! tangential DS, rotated by rad_c
         coef_r(2,3)=-sinD*sin2AAZ*rad_rot ! radial SS, rotated by rad_c
         coef_r(2,4)=cosAAZ*cosD*rad_rot ! radial DS, rotated by rad_c
         coef_r(2,5)=0 ! radial 45, rotated by rad_c

         return
         end

