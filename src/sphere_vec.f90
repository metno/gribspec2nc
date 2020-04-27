module sphere_vec

  use, intrinsic :: iso_fortran_env, only : real32, real64

  integer, parameter :: sng = real32
  integer, parameter :: dbl = real64

  real(kind=dbl), parameter :: pi     = acos(-1.0_dbl)
  real(kind=dbl), parameter :: halfpi = 0.5_dbl*pi
  real(kind=dbl), parameter :: twopi  = 2.0_dbl*pi
  real(kind=dbl), parameter :: rad    = pi/180.0_dbl
  real(kind=dbl), parameter :: rEarth = 6371000.0_dbl

contains

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
subroutine spheredist_vec(lon1, lat1, lon2, lat2, dist, nelem)

  implicit none

  real(kind=sng), intent(in)  :: lon1(nelem) ! [deg]
  real(kind=sng), intent(in)  :: lat1(nelem) ! [deg]
  real(kind=sng), intent(in)  :: lon2        ! [deg]
  real(kind=sng), intent(in)  :: lat2        ! [deg]
  real(kind=sng), intent(out) :: dist(nelem)
  integer,        intent(in)  :: nelem

  integer i
  real(kind=dbl) dTmp(nelem), lan2r, lat2r

  lan2r = real(lon2,kind=dbl)*rad
  lat2r = real(lat2,kind=dbl)*rad

  call spherearcrad_vec(real(lon1,kind=dbl)*rad, real(lat1,kind=dbl)*rad, lan2r, lat2r, dTmp, nelem)
  dist = real(dTmp*rEarth, kind=sng)

end subroutine spheredist_vec

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
subroutine spherearcrad_vec(rlon1, rlat1, rlon2, rlat2, dist, nelem)

  implicit none

  real(kind=dbl), intent(in)  :: rlon1(nelem) ! [rad]
  real(kind=dbl), intent(in)  :: rlat1(nelem) ! [rad]
  real(kind=dbl), intent(in)  :: rlon2        ! [rad]
  real(kind=dbl), intent(in)  :: rlat2        ! [rad]
  real(kind=dbl), intent(out) :: dist(nelem)  ! [rad]
  integer,        intent(in)  :: nelem

  integer i

  real(kind=dbl) a(nelem), d(nelem), l1(nelem), l2, l2c, l2s, s2
  real(kind=dbl) x1(nelem), y1(nelem), z1(nelem), x2, y2, z2 ! Cartesian positions

  l2  = angpi(rlon2)
  l2c = cos(l2)
  l2s = sin(l2)

  ! (x,y,z) of pos 2.
  s2 = sin(halfpi - rlat2)
  z2 = cos(halfpi - rlat2)
  x2 = s2*l2c
  y2 = s2*l2s

  do i=1,nelem
    l1(i) = angpi(rlon1(i))
  end do

  a(1:nelem)  = halfpi - rlat1(1:nelem)
  x1(1:nelem) = sin(a(1:nelem))*cos(l1(1:nelem)) ! (x,y,z) of pos 1.
  y1(1:nelem) = sin(a(1:nelem))*sin(l1(1:nelem))
  z1(1:nelem) = cos(a(1:nelem))

  d(1:nelem) = x1(1:nelem)*x2 + y1(1:nelem)*y2 + z1(1:nelem)*z2
  d(1:nelem) = min( 1.0_dbl, d(1:nelem))
  d(1:nelem) = max(-1.0_dbl, d(1:nelem))
  dist(1:nelem) = acos(d(1:nelem)) ! Arc length [rad]

end subroutine spherearcrad_vec

!!!!! FUNCTION angpi !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Maps arbitrary angle [rad] to [-PI, PI) [rad].
! Use this whenever two angles are subtracted.
!
! Requires ang2pi.
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=dbl) pure function angpi(ang)

  implicit none

  real(kind=dbl), intent(in) :: ang ! [rad]

  angpi = ang2pi(ang)
  angpi = angpi - pi*(sign(1.0_dbl, angpi-pi) + 1.0_dbl)
  ! angpi = modulo(ang+pi, twopi)-pi

end function angpi

!!!!! FUNCTION ang2pi !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Maps arbitrary angle [rad] to [0, 2*PI) [rad].
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=dbl) pure function ang2pi(ang)

  implicit none

  real(kind=dbl), intent(in) :: ang ! [rad]

  ang2pi = mod(ang, 2*pi) - (sign(1.0_dbl, ang) - 1.0_dbl)*pi

end function ang2pi

end module sphere_vec
