module sphere_vec

  use, intrinsic :: iso_fortran_env, only : real32, real64

  integer, parameter :: sng = real32
  integer, parameter :: dbl = real64

  real(kind=dbl), parameter :: pi     = acos(-1.0_dbl)
  real(kind=dbl), parameter :: halfpi = 0.5_dbl*pi
  real(kind=dbl), parameter :: rad    = pi/180.0_dbl
  real(kind=dbl), parameter :: rEarth = 6371000.0_dbl

  public :: spheredist, spherearcrad

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
real(kind=sng) function spheredist(lon1, lat1, lon2, lat2)

  implicit none

  real(kind=sng), intent(in) :: lon1 ! [deg]
  real(kind=sng), intent(in) :: lat1 ! [deg]
  real(kind=sng), intent(in) :: lon2 ! [deg]
  real(kind=sng), intent(in) :: lat2 ! [deg]

  spheredist = real(spherearcrad(lon1*rad,lat1*rad,lon2*rad,lat2*rad)*rEarth,kind=sng)

end function spheredist

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
real(kind=dbl) function spherearcrad(rlon1, rlat1, rlon2, rlat2)

  implicit none

  real(kind=dbl), intent(in) :: rlon1 ! [rad]
  real(kind=dbl), intent(in) :: rlat1 ! [rad]
  real(kind=dbl), intent(in) :: rlon2 ! [rad]
  real(kind=dbl), intent(in) :: rlat2 ! [rad]

  real(kind=dbl) a, c, d, l1, l2
  real(kind=dbl) x1, y1, z1, x2, y2, z2 ! Cartesian positions

  l1 = angpi(rlon1)
  a  = halfpi-rlat1

  l2 = angpi(rlon2)
  c  = halfpi-rlat2

  x1 = sin(a)*cos(l1) ! (x,y,z) of pos 1.
  y1 = sin(a)*sin(l1)
  z1 = cos(a)

  x2 = sin(c)*cos(l2) ! (x,y,z) of pos 2.
  y2 = sin(c)*sin(l2)
  z2 = cos(c)

  d = x1*x2 + y1*y2 + z1*z2
  d = min( 1.0_dbl, d)
  d = max(-1.0_dbl, d)

  spherearcrad = acos(d) ! Arc length [rad]

end function spherearcrad

!!!!! FUNCTION angpi !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Maps arbitrary angle [rad] to [-PI, PI) [rad].
! Use this whenever two angles are subtracted.
!
! Requires ang2pi.
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=dbl) function angpi(ang)

  implicit none

  real(kind=dbl), intent(in) :: ang ! [rad]

  angpi = ang2pi(ang)
  angpi = angpi - pi*(sign(1.0_dbl, angpi-pi) + 1.0_dbl)

end function angpi

!!!!! FUNCTION ang2pi !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Maps arbitrary angle [rad] to [0, 2*PI) [rad].
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=dbl) function ang2pi(ang)

  implicit none

  real(kind=dbl), intent(in) :: ang ! [rad]

  ang2pi = mod(ang, 2*pi) - (sign(1.0_dbl, ang) - 1.0_dbl)*pi

end function ang2pi

end module sphere_vec
