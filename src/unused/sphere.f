!!!!! CONTENTS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! REAL FUNCTION spheredist(lon1,lat1,lon2,lat2) 
! REAL FUNCTION spherearc(lon1,lat1,lon2,lat2)
! REAL FUNCTION spheredir(lon1,lat1,lon2,lat2) 
! SUBROUTINE    spherepos(lon1,lat1,distance,direction,lon2,lat2)
! REAL FUNCTION ang180(angle)
! REAL FUNCTION ang360(angle)
!
! Private double precision functions
!
! Computations are done in double precision, I/O in single precision
!
! 2001-08-23, Oyvind.Breivik@dnmi.no
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! FUNCTION spheredist !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Returns the great circle (geodesic) distance measured in meters
! between positions (lon1,lat1) and (lon2,lat2) [deg].
!
! Reference: Abramowitz and Stegun (1970): Handbook of mathematical functions,
! Section 4.3.149.
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
!
! See also spherearc.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION spheredist(lon1,lat1,lon2,lat2)

      IMPLICIT NONE

      REAL lon1, lat1, lon2, lat2  ! [deg]

      REAL spheredist              ! [deg]

      ! Functions called
      DOUBLE PRECISION spherearcrad

      ! Constants
      DOUBLE PRECISION RAD

      ! Locals
      REAL rearth0
      DOUBLE PRECISION rearth

      RAD=acos(-1.0D0)/180.0D0

      call earthr(rearth0)
      rearth = dble(rearth0)

      spheredist = real(spherearcrad(lon1*RAD,lat1*RAD,lon2*RAD,
     +                  lat2*RAD)*rearth)

      END ! FUNCTION spheredist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! FUNCTION spherearc !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Returns the great circle (geodesic) distance measured in degrees [0,180) 
! between positions (lon1,lat1) and (lon2,lat2) [deg].
!
! Reference: Abramowitz and Stegun (1970): Handbook of mathematical functions,
! Section 4.3.149.
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
!
! See also spheredist.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION spherearc(lon1,lat1,lon2,lat2)

      IMPLICIT NONE

      REAL lon1, lat1, lon2, lat2  ! [deg]

      REAL spherearc               ! [deg]

      ! Function called
      DOUBLE PRECISION spherearcrad

      ! Constants
      DOUBLE PRECISION DEG, RAD

      RAD=acos(-1.0D0)/180.0D0
      DEG=1.0D0/RAD

      spherearc = real(spherearcrad(lon1*RAD,lat1*RAD,lon2*RAD,lat2*RAD)
     +                 *DEG)

      END ! FUNCTION spherearc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! SUBROUTINE spherepos !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Finds the point (lon2,lat2) on the sphere separated from point (lon1,lat1) by ! dist [deg], a great circle path which crosses through points (lon1,lat1) and
! (lon2,lat2). This great circle path has local direction (heading) dir [deg]
! relative to north in position (lon1,lat1).
!
! Note: distance is in meters if dist<0.
!
! Reference: Abramowitz and Stegun (1970): Handbook of mathematical functions,
! Section 4.3.149.
!
! 2001-08-15, Oyvind.Breivik@dnmi.no.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE spherepos(lon1,lat1,dist0,dir,lon2,lat2)

      IMPLICIT NONE

      ! Interface
      REAL lon1, lat1, dist0, dir ! [deg]
      REAL lon2, lat2             ! [deg]

      ! Functions called
      REAL ang180, ang360

      ! Constants
      DOUBLE PRECISION PI, RAD, DEG, EPS

      ! Locals
      REAL rearth0
      DOUBLE PRECISION a, b, c, bb, cc, x, y, dist, rearth, ci

      ! Main

      PI=acos(-1.0D0)
      RAD=PI/180.0D0
      DEG=1.0D0/RAD
      EPS=1.0D-15

      call earthr(rearth0)
      rearth = dble(rearth0)
      if (dist0 .LT. 0.0) then
         dist=-dble(dist0)*DEG/rearth
      else
         dist=dble(dist0)
      endif
         
      a = RAD*(90.0D0-dble(lat1))
      b = RAD*dist
      cc = RAD*(360.0D0-dble(dir))

      ci = cos(a)*cos(b) + sin(a)*sin(b)*cos(cc)
      ci = min(1.0d0,ci)
      ci = max(-1.0d0,ci)
      c = acos(ci)
      y = sin(b)*sin(cc) / max(sin(c),EPS)
      x = (cos(b) - cos(c)*cos(a)) / max(sin(a)*sin(c),EPS)
      bb = atan2(y,x)

      lon2 = ang180(lon1-real(DEG*bb))
      lat2 = ang180(90.0-real(DEG*c))

      END ! SUBROUTINE spherepos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! FUNCTION spheredir !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Returns the local direction measured clockwise relative to north from
! point (lon1,lat1) towards (lon2,lat2) following a great circle
! (geodesic) path on a sphere.
!
! Requires spherearcrad, ang180, ang360.
!
! Reference: Abramowitz and Stegun (1970): Handbook of mathematical functions,
! Section 4.3.149.
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION spheredir(lon1,lat1,lon2,lat2)

      IMPLICIT NONE

      REAL lon1, lat1, lon2, lat2 ! [deg]

      REAL spheredir              ! [deg]

      ! Functions called
      DOUBLE PRECISION spherearcrad
      REAL ang180, ang360

      ! Constants
      DOUBLE PRECISION RAD, DEG, EPS

      ! Locals
      DOUBLE PRECISION a, b, c, bb, cc, x, y

      RAD=acos(-1.0D0)/180.0D0
      DEG=1.0D0/RAD
      EPS=1.0D-15

      a = RAD*(90.0D0-lat1)
      c = RAD*(90.0D0-lat2)
      bb = RAD*ang180(lon1-lon2)
   
      b = spherearcrad(lon1*RAD,lat1*RAD,lon2*RAD,lat2*RAD)
   
      y = sin(c)*sin(bb)/sin(b)
      x = (cos(c) - cos(a)*cos(b)) / max(sin(a)*sin(b),EPS)

      cc = atan2(y,x)

      spheredir = ang360(360.0-real(DEG*cc))

      END ! FUNCTION spheredir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! FUNCTION ang180 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Maps arbitrary angle [deg] to [-180, 180) [deg].
! Use this whenever two angles are subtracted.
! Requires ang360.
!
! 2001-08-15, Oyvind.Breivik@dnmi.no.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION ang180(ang)

      IMPLICIT NONE

      REAL ang    ! [deg]

      REAL ang180 ! [deg]

      REAL ang360 ! [deg]

      ang180 = ang360(ang)
      ang180 = ang180 - 180.0*(sign(1.0,ang180-180.0)+1.0)

      END ! FUNCTION ang180
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! FUNCTION ang360 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Maps arbitrary angle [deg] to [0, 360) [deg].
!
! 2001-08-15, Oyvind.Breivik@dnmi.no.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION ang360(ang)

      IMPLICIT NONE

      REAL ang    ! [deg]

      REAL ang360 ! [deg]

      ang360 = mod(ang, 360.0) - (sign(1.0,ang)-1.0)*180.0

      END ! FUNCTION ang360
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! FUNCTION spherearcrad !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Returns the great circle (geodesic) distance measured in radians [0,PI)
! between positions (rlon1,rlat1) and (rlon2,rlat2) [rad].
!
! Requires angpi.
!
! Reference: Abramowitz and Stegun (1970): Handbook of mathematical functions,
! Section 4.3.149.
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
!
! Note the difference from spherdist which returns the distance in meters.
!
! See also spherearc, spherearc2 and spherdist.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION spherearcrad(rlon1,rlat1,rlon2,rlat2)

      IMPLICIT NONE

      DOUBLE PRECISION rlon1, rlat1, rlon2, rlat2  ! [rad]

      DOUBLE PRECISION spherearcrad                ! [rad]

      !!! Functions
      DOUBLE PRECISION angpi

      !!! Locals

      DOUBLE PRECISION PI

      DOUBLE PRECISION a, c, d, l1, l2
      DOUBLE PRECISION x1, y1, z1, x2, y2, z2          ! Cartesian positions

      PI=acos(-1.0D0)

      l1 = angpi(rlon1)
      a  = 0.5D0*PI-rlat1

      l2 = angpi(rlon2)
      c  = 0.5D0*PI-rlat2 

      x1 = sin(a)*cos(l1)                     ! (x,y,z) of pos 1.
      y1 = sin(a)*sin(l1)
      z1 = cos(a)

      x2 = sin(c)*cos(l2)                     ! (x,y,z) of pos 2.
      y2 = sin(c)*sin(l2)
      z2 = cos(c)
   
      d = x1*x2+y1*y2+z1*z2
      d = min(1.0d0,d)
      d = max(-1.0d0,d)
      spherearcrad = acos(d)  ! Arc length [rad]

      END ! FUNCTION spherearcrad
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! FUNCTION ang2pi !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Maps arbitrary angle [rad] to [0, 2*PI) [rad].
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION ang2pi(ang)

      IMPLICIT NONE

      DOUBLE PRECISION ang    ! [rad]

      DOUBLE PRECISION ang2pi ! [rad]

      DOUBLE PRECISION PI

      PI=acos(-1.0D0)

      ang2pi = mod(ang, 2*PI) - (sign(1.0D0,ang)-1.0D0)*PI

      END ! FUNCTION ang2pi 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!! FUNCTION angpi !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Maps arbitrary angle [rad] to [-PI, PI) [rad].  
! Use this whenever two angles are subtracted.  
!
! Requires ang2pi.
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION angpi(ang) 

      IMPLICIT NONE

      DOUBLE PRECISION ang   ! [rad]

      DOUBLE PRECISION angpi ! [rad]

      DOUBLE PRECISION ang2pi

      DOUBLE PRECISION PI

      PI=acos(-1.0D0)

      angpi = ang2pi(ang)
      angpi = angpi - PI*(sign(1.0D0,angpi-PI)+1.0D0)

      END ! FUNCTION angpi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!! SUBROUTINE earthr !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  NAME:
!     earthr
!
!  PURPOSE:
!     Return earth radius [m], the DNMI standard value
!
!  SYNOPSIS:
!     subroutine earthr(rearth)
!     real rearth
!
!  OUTPUT:
!     rearth  - earth radius i unit meter
!
!-----------------------------------------------------------------------
!  DNMI/FoU  xx.xx.19xx  ............ ... rearth = 6368.00 km in models
!  DNMI/FoU  25.08.1995  Anstein Foss ... rearth = 6371.22 km
!  DNMI/FoU  21.03.1996  Anstein Foss ... rearth = 6371.00 km
!-----------------------------------------------------------------------
!
      SUBROUTINE earthr(rearth)

      IMPLICIT NONE

      REAL rearth ! [m]

      rearth = 6371000.0

      return
      END ! SUBROUTINE rearth
