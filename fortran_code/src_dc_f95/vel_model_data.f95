module vel_model_data


   use constants, only : nlay, epsilon
   implicit none
   integer, private :: n_layers
   integer :: n_layers_new, src, rcv
   real, private :: vel_p(nlay), vel_s(nlay), dens(nlay), thick(nlay), qp(nlay), qs(nlay)
   real*8 :: new_vel_p(nlay), new_vel_s(nlay), new_dens(nlay), new_thick(nlay), qqp(nlay), qqs(nlay)
   real*8 :: slowness(nlay), xi(nlay), mu(nlay), depths


contains


   subroutine read_vel_model(vel_model)
!
!  Args:
!  vel_model: file with velocity model
!
   implicit none
   character(len=100) :: vel_model
   integer :: j
   vel_model = trim(vel_model)
   open(1, file=vel_model, status='old')
   read(1, *) n_layers
   do j = 1, n_layers
      read (1, *)vel_p(j), vel_s(j), dens(j), thick(j), qp(j), qs(j)
   enddo
   if (thick(n_layers).gt.1.e-10) then
      n_layers = n_layers + 1
      thick(n_layers) = 0.0
      vel_p(n_layers) = vel_p(n_layers - 1)
      vel_s(n_layers) = vel_s(n_layers - 1)
      dens(n_layers) = dens(n_layers - 1)
      qp(n_layers) = qp(n_layers - 1)
      qs(n_layers) = qs(n_layers - 1)
   endif
   end subroutine read_vel_model

  
   subroutine update_model(depth_source, depth_rcv)
!
!  Args:
!  depth_source: depth of point source
!  depth_rcv: depth of receiver
!
   implicit none
   real, intent(in) :: depth_source, depth_rcv
   real :: cum_dep, depths(nlay), depth1, depth2
   integer :: j
   depth1 = min(depth_source, depth_rcv)
   depth2 = max(depth_source, depth_rcv)
   depths(1) = 0.0
   cum_dep = 0.0
   new_vel_p = vel_p
   new_vel_s = vel_s
   new_dens = dens
   new_thick = thick
   qqp = qp
   qqs = qs
   do j = 1, n_layers
      cum_dep = cum_dep + thick(j)
      if (thick(j).gt.1.e-6) then
         depths(j + 1) = cum_dep
      else
         depths(j + 1) = 99999.
      endif
   enddo
   n_layers_new = n_layers

   rcv = 1
   do j = 1, n_layers
      if((depth1 .gt. depths(j)) .and. (depth1 .le. depths(j + 1))) then
         rcv = j + 1
         exit
      end if
   enddo

   if (rcv.gt.1) then
      new_thick(rcv - 1) = depth1 - depths(rcv - 1)
      new_thick(rcv:n_layers + 1) = thick(rcv - 1:n_layers)
      new_vel_p(rcv:n_layers + 1) = vel_p(rcv - 1:n_layers)
      new_vel_s(rcv:n_layers + 1) = vel_s(rcv - 1:n_layers)
      new_dens(rcv:n_layers + 1) = dens(rcv - 1:n_layers)
      qqp(rcv:n_layers + 1) = qp(rcv - 1:n_layers)
      qqs(rcv:n_layers + 1) = qs(rcv - 1:n_layers)
      new_thick(rcv) = depths(rcv) - depth1
      if(new_thick(rcv) .lt. 1.e-6)new_thick(rcv) = 1.e-6
      if(new_thick(rcv) .ge. 6371.)new_thick(rcv) = 0.0
      n_layers_new = n_layers_new + 1
   endif
 
   src = 1
   do j = 1, n_layers
      if((depth2 .gt. depths(j)) .and. (depth2 .le. depths(j + 1))) then
         src = j + 1
         exit
      end if
   enddo

   if (src.gt.1) then
      new_thick(src - 1) = depth2 - depths(src - 1)
      new_thick(src:n_layers + 1) = thick(src - 1:n_layers)
      new_vel_p(src:n_layers + 1) = vel_p(src - 1:n_layers)
      new_vel_s(src:n_layers + 1) = vel_s(src - 1:n_layers)
      new_dens(src:n_layers + 1) = dens(src - 1:n_layers)
      qqp(src:n_layers + 1) = qp(src - 1:n_layers)
      qqs(src:n_layers + 1) = qs(src - 1:n_layers)
      new_thick(src) = depths(src) - depth2
      if(new_thick(src) .lt. 1.e-6)new_thick(src) = 1.e-6
      if(new_thick(src) .ge. 6371.)new_thick(src) = 0.0
      n_layers_new = n_layers_new + 1
   endif

   if (depth_source.lt.depth_rcv) then
      j = rcv
      rcv = src
      src = j
   endif

   do j = 1, n_layers_new
      if (new_vel_p(j).lt.0.001) then
         slowness(j) = 0.001**2
      else
         slowness(j) = new_vel_p(j)**2
      endif
   enddo
   slowness = 1./slowness
!   slowness = 1. / new_vel_p ** 2

   end subroutine update_model


   subroutine flip_model()
   implicit none
   real*8 :: temp(nlay, 6)
   integer :: i, j
   
   do i=1, n_layers_new
      j = n_layers_new-i+1
      temp(i,1)=new_vel_p(i)
      temp(i,2)=new_vel_s(i)
      temp(i,3)=new_dens(i)
      temp(i,4)=qqp(i)
      temp(i,5)=qqs(i)
      temp(i,6)=new_thick(i)
   enddo
   new_vel_p(:)=temp(:,1)
   new_vel_s(:)=temp(:,2)
   new_dens(:)=temp(:,3)
   qqp(:)=temp(:,4)
   qqs(:)=temp(:,5)
   new_thick(:)=temp(:,6)
    
   end subroutine flip_model


   subroutine extra(hs)
!
!  Args:
!  hs: depth of source
!
   implicit none
   real*8 :: hs
   integer :: i

   do i=1,n_layers_new
      if (new_vel_s(i).lt.epsilon) new_vel_s(i)=epsilon
      if ( i.lt.src .and. i.ge.rcv ) hs = hs + new_thick(i)
   enddo
   xi = new_vel_s*new_vel_s/(new_vel_p*new_vel_p)
   mu = new_dens*new_vel_s*new_vel_s
   end subroutine extra
   
   
   subroutine extra2(hs)
!
!  Args:
!  hs: depth of source
!
   implicit none
   real*8 :: hs
   integer :: i

   do i=1,n_layers_new
      write(0,*) 'Input thickness Vp Vs rho Qa Qb for layer',i
      write(0,'(i4,4f7.2,2e9.2)')i,new_thick(i),new_vel_p(i),new_vel_s(i),new_dens(i),qqp(i),qqs(i)
   enddo
   write(0,'(a15,f8.3)')'source-station separation=',hs

   end subroutine extra2


end module vel_model_data
