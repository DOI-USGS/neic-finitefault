!       Make the Green Function bank 
!       Input dist_max dist_min d_step
!       Input dep_max dep_min dep_step
!       Input Lnpt dt
!       Input Green_name
!       output Green
program green_bank_fk


   use constants, only : nlay, ndis, nt
   use retrieve_gf, only : block_gg, dt, lnpt, dist_max, dist_min, d_step, dep_max, dep_min, dep_step, get_gf_data, t_cor
   use vel_model_data, only : read_vel_model, update_model, src, rcv 
   use wave_travel, only : trav_fk
   use bessel2, only : load_bessel
   use fk, only : sub_bs_dc
   implicit none
   real dist(ndis)
   real t0(ndis),tmin(ndis),depth
   real :: green(nt, 8, ndis)
   integer iz, k, ll, n_com, nd_max, npt, ntc, nx, nz
   logical :: disp

   character(len=100) gf_file, vel_model, gf_bank, input
!
   call getarg(1, input)
   disp = (input.eq.'cgnss')

   gf_file = 'Green.in'
   if (disp) gf_file = 'Green_cgnss.in'
   call get_gf_data(gf_file, vel_model, gf_bank)
   npt=2**lnpt
!
   if(abs(dist_min-int(dist_min/d_step+0.001)*d_step).gt.1.0e-4)then
      write(*,*) "To improve the speed of this calculation, "
      write(*,*) "The minimum distance dist_min = m *d_step "
      write(*,*)dist_min,d_step,abs(dist_min-int(dist_min/d_step+0.001)*d_step)
!    pause
   endif
   if(abs(dist_max-int(dist_max/d_step+0.001)*d_step).gt.1.0e-4)then
      write(*,*) "To improve the speed of this calculation, "
      write(*,*) "The maximum distance dist_min = m *d_step "
      write(*,*)dist_max,d_step,abs(dist_max-int(dist_max/d_step+0.001)*d_step)
!	pause
   endif
   nd_max=int(dist_max/d_step+0.1)+1

   if(npt.gt.nt)then
      write(0,*)"The length of green's function must be shorten"
      stop
   endif
   if(dep_max.lt.dep_min)then
      write(*,*)"The depth region is wrong"
      stop
   endif
   if(dist_max.lt.dist_min)then
      write(*,*)"distance region is wrong"
      stop
   endif

   nx=int((dist_max-dist_min)/d_step)+1
   nz=int((dep_max-dep_min)/dep_step)+1

   write(*,*)'total =',nx,nz
!	open(100,file='stored_bessel')
!	open(101,file='computed_bessel')
   if(nx.gt.ndis) then
      write(*,*)"please reduce the number of Green function"
      stop
   endif
   open(11,file=gf_bank,status='unknown',access='direct',recl=block_gg)

   call read_vel_model(vel_model)

   do k=1,nx
      dist(k)=dist_min+(k-1)*d_step
   enddo
   call load_bessel(dist, nd_max)
   ll=0
   do iz=1,nz
      depth=dep_min + (iz-1)*dep_step
      write(*,*) depth
      call update_model(depth,0.0)
      call trav_fk(dist,tmin,nx)
      write(*,*)'t0: ', (tmin(k),k=1,20)
      do k=1,nx
         t0(k)=tmin(k)-t_cor
      enddo
      write(0,*)src, rcv
      call sub_bs_dc(nx,dist,t0,green,disp)
      ll = nx * (iz - 1)
      do k=1,nx
         ll=ll+1
         write(11,rec=ll) &
            iz,k,dist(k),t0(k),depth,dt,npt,((green(ntc,n_com,k),ntc=1,npt),n_com=1,8)
      enddo
   enddo
   close(11)

 
end program green_bank_fk

