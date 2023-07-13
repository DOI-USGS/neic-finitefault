module geodesics


   use constants, only : pi
   implicit none
   real*8, private :: c0, c2, c4
   real*8, parameter, private :: rad = 6378.155d0, const = 57.2957795d0, gl = .99327733d0, ec = 0.672267002d-2
   real*8, private :: th(2), phi(2), xdeg(2), azinv(2), dist(2), az(2)
   public  :: distaz
   private :: disaz1, arc_dist


contains


   subroutine distaz(lat_sta, lon_sta, lat_e, lon_e, dis, azz, baz)
!
!  Args:
!  lat_sta: station latitude
!  lon_sta: station longitude
!  lat_e: point_source latitude
!  lon_e: point source longitude
!  dis: source-station distance
!  azz: source-station azimuth
!  baz: source-station back-azimuth
!
   implicit none
   real, intent(in) :: lat_sta, lon_sta, lat_e, lon_e
   real*8, intent(out) :: dis, azz, baz
   if(abs(lat_sta-lat_e).gt.1.0e-5 .or. abs(lon_sta-lon_e).gt.1.e-5)then
      th(1) = lat_sta
      th(2) = lat_e
      phi(1) = lon_sta
      phi(2) = lon_e
      call disaz1(1, 2)
      azz = az(1)
      baz = azinv(1)
      dis = dist(1)
   else 
      dis = 0.d0
      azz = 0.d0
      baz = 180.d0
   endif
   end subroutine distaz


   subroutine disaz1(n, k)
!
!  Args:
!  n: amount of stations for which to compute distance, azimuth and back-azimuth
!  k: index of point source
!
   implicit none
   integer, intent(in) :: n, k
   integer :: i
   real*8 :: sink, cosk, tank, sini, cosi, tani
   real*8 :: a, b, c, a1, b1, c1, sc, sd1, sd, arg
   real*8 :: el, x, cdl, sdl, dum, a12, dum1, cs12
   real*8 :: v1, v2, z1, z2, x2, y2, chck, u1, u2, b0
   real*8 :: test1, test2, test3, coeff1, coeff2, coeff3
   real*8 :: al12, al21, eo

   th(:) = th(:) / const
   phi(:) = phi(:) / const   
   sink = sin(th(k))
   cosk = cos(th(k))
   tank = sink / cosk
   arg = atan(gl * tank)
   c = sin(arg)
   a = cos(arg) * cos(phi(k))
   b = cos(arg) * sin(phi(k))
   do i = 1, n
      sini = sin(th(i))
      cosi = cos(th(i))
      tani = sini / cosi
      arg = atan(gl * tani)
      C1 = sin(arg)
      A1 = cos(arg) * cos(phi(i))
      B1 = cos(arg) * sin(phi(i))
      sc =  a * a1 + b * b1 + c * c1
      sd1 = (a - a1) ** 2 + (b - b1) ** 2 + (c - c1) ** 2
      sd = (a + a1) ** 2 + (b + b1) ** 2 + (c + c1) ** 2
      sd = sqrt(sd1 * sd  / 4.d0)
      xdeg(i) = atan(sd / sc) * const
      if (sc.lt.0) xdeg(i) = xdeg(i) + 180.d0
      el = ec / (1.d0 - ec)
      x = (1.d0 + el + tani ** 2) / (1.d0 + el + TANK ** 2)
      CDL = cos(phi(i) - phi(k))
      SDL = sin(phi(i) - phi(k))
      al12 = tani / (1.d0 + el) / tank + ec * sqrt(x) - cdl
      al21 = tank / (1.d0 + el) / tani + ec * sqrt(1.d0 / x) - cdl
      DUM = sink * al12
      A12 = atan(sdl / dum)
      DUM1 = sini * al21
      az(i) = A12 * const
      azinv(i) = atan(- sdl / dum1) * const
      if (sdl.lt.0) then
         if (dum .lt. 0)az(i) = az(i) - 180.d0
      else
         if (dum .lt. 0)az(i) = 180.d0 + az(i)
      endif
      if (-sdl.lt.0) then
         if (dum1.lt.0)azinv(i) = azinv(i) - 180.d0
      else
         if (dum1.lt.0)azinv(i) = azinv(i) + 180.d0
      endif
      if(az(i).lt.0.0) az(i) = 360.d0 + az(i)
      if(azinv(i).lt.0.0) azinv(i) = 360.d0 + azinv(i)
      cs12 = cos(a12)
      th(i) = th(i) * const
      phi(i) = phi(i) * const
      eo = EL * ((COSK * CS12) ** 2 + SINK ** 2)
      c0 = 1.d0 + eo / 4.0 - 3.d0 * eo ** 2 / 64.d0 + 5.d0 * eo ** 3 / 256.d0
      c2 = - eo / 8.d0 + eo ** 2 / 32.d0 -15.d0 * eo ** 3 / 1024.d0
      c4 = - eo ** 2 / 256.d0 + 3.d0 * eo ** 3 / 1024.d0
      V1 = rad / sqrt(1.d0 - ec * sink ** 2)
      V2 = rad / sqrt(1.d0 - ec * sini ** 2)
      Z1 = V1 * (1.d0 - EC) * SINK
      Z2 = V2 * (1.d0 - EC) * SINI
      X2 = V2 * COSI * CDL
      Y2 = V2 * COSI * SDL

      chck = xdeg(i) * 111.32d0
      arg = TANK / CS12 / sqrt(1.d0 + EO)
      u1 = atan(arg)
      X = (X2 * CS12 - Y2 * SINK * sin(A12)) * sqrt(1.d0 + EO)
      arg = (V1 * SINK + (1.d0 + EO) * (Z2 - Z1)) / X
      u2 = atan(arg)
      x = 1.d0 + EL * (COSK * CS12) ** 2
      b0 = V1 * sqrt(x) / (1.d0 + EO)
      coeff1 = arc_dist(b0, u1, u2)
      test1 = abs(chck - abs(coeff1))
      coeff2 = arc_dist(b0, u1, u2 + pi)
      test2 = abs(chck - abs(coeff2))
      coeff3 = arc_dist(b0, u1, u2 - pi)
      test3 = abs(chck - abs(coeff3))
      if(test1.lt.30) then
         dist(i) = abs(coeff1)
      elseif(test2.lt.30) then
         dist(i) = abs(coeff2)
      elseif(test3.lt.30) then
         dist(i) = abs(coeff3)
      else
         dist(i) = chck
      endif
   enddo
   end subroutine disaz1


   function arc_dist(b, u, s) result(dist)
!
!  Args:
!  b: ?
!  u: ?
!  s: ?
!
   implicit none
   real*8 :: b, u, s, dist, x, sinu2, sinu4, sins2, sins4
   sinu2 = sin(2.d0 * u)
   sinu4 = sin(4.d0 * u)
   sins2 = sin(2.d0 * s)
   sins4 = sin(4.d0 * s)
   x = c0 * (u - s) + c2 * (sinu2 - sins2) + c4 * (sinu4 - sins4)
   dist = b * x
   end function arc_dist


end module geodesics
