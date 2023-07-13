!
!	lupei's format. the source and receiver at the top of correponding layers
!
module cmodel


contains


   subroutine cmodel_fk(depth_s,depth_r, jjo, nlay_s, nlay_r, c, s, den, th, qa, qb)
!
!  Args:
!  depth_s: depth of source
!  depth_r: depth of receiver
!  jjo: layers of modified velocity model
!  nlay_s: layer of source
!  nlay_r: layer of receiver
!  c: P wave velocity of modified velocity model
!  s: S wave velocity of modified velocity model
!  den: density of modified velocity model
!  th: thickness of each layer for modified velocity model
!  qa: attenuation for modified velocity model
!  qb: attenuation for modified velocity model
!
   implicit none
   include 'fk.h'
   integer :: io_up
   integer :: jo,jjo,nb,j,nlay_s,nlay_r
   real :: depth_s,depth_r,depth_1,depth_2
   real :: c(0:nlay), s(0:nlay),den(0:nlay),th(0:nlay),qa(0:nlay),qb(0:nlay)
   real :: x,cc(0:nlay),ss(0:nlay),dd(0:nlay),tth(0:nlay),dep(0:nlay)
   real :: qqa(0:nlay),qqb(0:nlay)
   common/model_fk/jo,cc,ss,dd,tth, qqa, qqb 
!
!	Read in velocity model and split source layer
!
   depth_1=depth_r
   depth_2=depth_s
   io_up=0
   if(depth_s.lt.depth_r)then
      depth_1=depth_s
      depth_2=depth_r
      io_up=1
   endif
   x=0.0
   dep(1)=0.0
   jjo = jo
   do j=1,jjo
      x=x+tth(j)
      if(tth(j).gt.1e-6)then
         dep(j+1)=x
      else
!	  half space
         dep(j+1)=99999.
      endif
      c(j)=cc(j)
      s(j)=ss(j)
      den(j)=dd(j)
      th(j)=tth(j)
      qa(j)=qqa(j)
      qb(j)=qqb(j)
   enddo
   nb=1
   do j=1,jjo
      if(depth_1.gt.dep(j).and.depth_1.le.dep(j+1)) then
         nb=j+1
         exit
      end if
   enddo
   if(nb.gt.1)then
      th(nb-1)=depth_1-dep(nb-1)
      do j=jjo+1,nb,-1
         c(j)=c(j-1)
         s(j)=s(j-1)
         den(j)=den(j-1)
         th(j)=th(j-1)
         qa(j)=qa(j-1)
         qb(j)=qb(j-1)
      enddo
      th(nb)=dep(nb)-depth_1
      if(th(nb).eq.0.0)th(nb)=0.001
      if(th(nb).gt.6371.0)th(nb)=0.0
      jjo=jjo+1
   endif
   if(io_up.eq.0)then
      nlay_r=nb
   else
      nlay_s=nb
   endif
   x=0.0
   dep(1)=0.0
   do j=1,jjo
      x=x+th(j)
      if(th(j).gt.1e-6)then
         dep(j+1)=x
      else
!        half space
         dep(j+1)=99999.
      endif
   enddo
   nb=1
   do j=1,jjo
      if(depth_2.gt.dep(j).and.depth_2.le.dep(j+1)) then
         nb=j+1
         exit
      end if
   enddo
   if(nb.gt.1)then
      th(nb-1)=depth_2-dep(nb-1)
      do j=jjo+1,nb,-1
         c(j)=c(j-1)
         s(j)=s(j-1)
         den(j)=den(j-1)
         th(j)=th(j-1)
         qa(j)=qa(j-1)
         qb(j)=qb(j-1)
      enddo
      th(nb)=dep(nb)-depth_2
      if(th(nb).eq.0.0)th(nb)=0.001
      if(th(nb).gt.6371.0)th(nb)=0.0
      jjo=jjo+1
   endif
   if(io_up.eq.0)then
      nlay_s=nb
   else
      nlay_r=nb
   endif
!       write(*,*)'end of cmodel'
   end subroutine cmodel_fk


end module cmodel

