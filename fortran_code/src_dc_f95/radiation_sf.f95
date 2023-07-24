!
! This program is modified from Lupei Zhu's code
!
!/******************************************************
!horizontal radiation coefficients of a single-force
!   In:
!        stk: strike_of_obs w.r.t to strike of the force
!                 measured clockwise (az_obs - az_sr)
!        dip: dip of the force, from horizontal down
!
!   algorithm:
!        vertical (UP) = f3*Z0 + (f1*cos(theta)+f2*sin(theta))*Z1
!        radial  (OUT) = f3*R0 + (f1*cos(theta)+f2*sin(theta))*R1
!        tangen   (CW) =       - (f1*sin(theta)-f2*cos(theta))*T1
!    where F = (0,cos(dip),-sin(dip))
!******************************************************/
!void    sf_radiat(float stk,float dip,float rad[4][3]) {
!   float sstk,sdip,cstk,cdip;
!   stk*=DEG2RAD; dip*=DEG2RAD;
!   sstk=sin(stk);cstk=cos(stk);
!   sdip=sin(dip);cdip=cos(dip);
!   rad[0][0]=-sdip;
!   rad[0][1]=rad[0][0];
!   rad[0][2]=0.;
!   rad[1][0]=cdip*cstk;
!   rad[1][1]=rad[1][0];
!   rad[1][2]=cdip*sstk;
!}
module radiation


contains


   subroutine sf_radiat(az,stk,dip,a)
!
!  Args:
!  az: source-station azimuth 
!  stk: strike of point source
!  dip: dip of point source
!  a: radiation pattern
!
   real stk, dip, az,a(3),const
   real sstk, ddip
   parameter(const=3.14159265358979/180.0)
   sstk=(az-stk)*const
   ddip=ddip*const
   a(1)=-sin(ddip)
   a(2)=cos(ddip)*sin(sstk)
   a(3)=cos(ddip)*cos(sstk)
   end subroutine sf_radiat


end module radiation      
