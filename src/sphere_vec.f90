
! ============================================================================ !
!
! Computations are done in double precision, I/O in single precision.
!
! 2001-08-23, Oyvind.Breivik@dnmi.no
! 2020-08-28, Veronica.Berglyd.Olsen@met.no :
!             Added vectorised versions and improved performance
!
! ============================================================================ !

module sphere_vec

  use, intrinsic :: iso_fortran_env, only : real32, real64

  implicit none

  integer, parameter :: sng = real32
  integer, parameter :: dbl = real64

  real(kind=dbl), parameter :: pi     = acos(-1.0_dbl)
  real(kind=dbl), parameter :: halfpi = 0.5_dbl*pi
  real(kind=dbl), parameter :: twopi  = 2.0_dbl*pi
  real(kind=dbl), parameter :: rad    = pi/180.0_dbl
  real(kind=dbl), parameter :: rEarth = 6371000.0_dbl

  interface spheredist
    module procedure spheredist_scalar
    module procedure spheredist_vector
  end interface spheredist

  interface spherearcrad
    module procedure spherearcrad_scalar
    module procedure spherearcrad_vector
  end interface spherearcrad

  private spheredist_scalar
  private spheredist_vector
  private spherearcrad_scalar
  private spherearcrad_vector

contains

! === SUBROUTINE spheredist ================================================== !
! Returns the great circle (geodesic) distance measured in meters
! between positions (lon1,lat1) and (lon2,lat2) [deg].
!
! Reference: Abramowitz and Stegun (1970): Handbook of mathematical functions,
! Section 4.3.149.
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
! 2020-04-28, Veronica.Berglyd.Olsen@met.no : Optimised and vectorised
!
! See also spherearc.
! ============================================================================ !
subroutine spheredist_scalar(lon1, lat1, lon2, lat2, dist)

  real(kind=sng), intent(in)  :: lon1 ! [deg]
  real(kind=sng), intent(in)  :: lat1 ! [deg]
  real(kind=sng), intent(in)  :: lon2 ! [deg]
  real(kind=sng), intent(in)  :: lat2 ! [deg]
  real(kind=sng), intent(out) :: dist

  real(kind=dbl) dTmp, lan2r, lat2r

  lan2r = real(lon2,kind=dbl)*rad
  lat2r = real(lat2,kind=dbl)*rad

  call spherearcrad_scalar(real(lon1,kind=dbl)*rad, real(lat1,kind=dbl)*rad, lan2r, lat2r, dTmp)
  dist = real(dTmp*rEarth, kind=sng)

end subroutine spheredist_scalar

subroutine spheredist_vector(lon1, lat1, lon2, lat2, dist, nelem)

  real(kind=sng), intent(in)  :: lon1(nelem) ! [deg]
  real(kind=sng), intent(in)  :: lat1(nelem) ! [deg]
  real(kind=sng), intent(in)  :: lon2        ! [deg]
  real(kind=sng), intent(in)  :: lat2        ! [deg]
  real(kind=sng), intent(out) :: dist(nelem)
  integer,        intent(in)  :: nelem

  real(kind=dbl) dTmp(nelem), lan2r, lat2r

  lan2r = real(lon2,kind=dbl)*rad
  lat2r = real(lat2,kind=dbl)*rad

  call spherearcrad_vector(real(lon1,kind=dbl)*rad, real(lat1,kind=dbl)*rad, lan2r, lat2r, dTmp, nelem)
  dist = real(dTmp*rEarth, kind=sng)

end subroutine spheredist_vector

! === SUBROUTINE spherearcrad ================================================ !
! Returns the great circle (geodesic) distance measured in radians [0,PI)
! between positions (rlon1,rlat1) and (rlon2,rlat2) [rad].
!
! Requires angpi.
!
! Reference: Abramowitz and Stegun (1970): Handbook of mathematical functions,
! Section 4.3.149.
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
! 2020-04-28, Veronica.Berglyd.Olsen@met.no : Optimised and vectorised
!
! Note the difference from spherdist which returns the distance in meters.
!
! See also spherearc, spherearc2 and spherdist.
! ============================================================================ !
subroutine spherearcrad_scalar(rlon1, rlat1, rlon2, rlat2, dist)

  real(kind=dbl), intent(in)  :: rlon1 ! [rad]
  real(kind=dbl), intent(in)  :: rlat1 ! [rad]
  real(kind=dbl), intent(in)  :: rlon2 ! [rad]
  real(kind=dbl), intent(in)  :: rlat2 ! [rad]
  real(kind=dbl), intent(out) :: dist

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

  dist = acos(d) ! Arc length [rad]

end subroutine spherearcrad_scalar

subroutine spherearcrad_vector(rlon1, rlat1, rlon2, rlat2, dist, nelem)

  real(kind=dbl), intent(in)  :: rlon1(nelem) ! [rad]
  real(kind=dbl), intent(in)  :: rlat1(nelem) ! [rad]
  real(kind=dbl), intent(in)  :: rlon2        ! [rad]
  real(kind=dbl), intent(in)  :: rlat2        ! [rad]
  real(kind=dbl), intent(out) :: dist(nelem)  ! [rad]
  integer,        intent(in)  :: nelem

  integer i

  real(kind=dbl) as(nelem), d(nelem), s2
  real(kind=dbl) x1(nelem), y1(nelem), z1(nelem), x2, y2, z2 ! Cartesian positions

  ! (x,y,z) of pos 2.
  s2 = sin(halfpi - rlat2)
  z2 = cos(halfpi - rlat2)
  x2 = s2*cos(rlon2)
  y2 = s2*sin(rlon2)

  ! (x,y,z) of pos 1.
  as(1:nelem) = sin(halfpi - rlat1(1:nelem))
  z1(1:nelem) = cos(halfpi - rlat1(1:nelem))
  x1(1:nelem) = as(1:nelem)*cos(rlon1(1:nelem))
  y1(1:nelem) = as(1:nelem)*sin(rlon1(1:nelem))

  d(1:nelem) = x1(1:nelem)*x2 + y1(1:nelem)*y2 + z1(1:nelem)*z2
  d(1:nelem) = min( 1.0_dbl, d(1:nelem))
  d(1:nelem) = max(-1.0_dbl, d(1:nelem))
  dist(1:nelem) = acos(d(1:nelem)) ! Arc length [rad]

end subroutine spherearcrad_vector

! === FUNCTION angpi ========================================================= !
! Maps arbitrary angle [rad] to [-PI, PI) [rad].
! Use this whenever two angles are subtracted.
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
! 2020-04-28, Veronica.Berglyd.Olsen@met.no : Rewritten to use modulo()
! ============================================================================ !
real(kind=dbl) pure function angpi(ang)
  real(kind=dbl), intent(in) :: ang ! [rad]
  angpi = modulo(ang + pi, twopi) - pi
end function angpi

! === FUNCTION ang2pi ======================================================== !
! Maps arbitrary angle [rad] to [0, 2*PI) [rad].
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
! 2020-04-28, Veronica.Berglyd.Olsen@met.no : Rewritten to use modulo()
! ============================================================================ !
real(kind=dbl) pure function ang2pi(ang)
  real(kind=dbl), intent(in) :: ang ! [rad]
  ang2pi = modulo(ang, twopi)
end function ang2pi

! === FUNCTION ang180 ======================================================== !
! Maps arbitrary angle [deg] to [-180, 180) [deg].
! Use this whenever two angles are subtracted.
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
! 2021-06-08, Veronica.Berglyd.Olsen@met.no : Moved from sphere.f
! ============================================================================ !
real pure function ang180(ang)

  implicit none
  real, intent(in) :: ang ! [deg]

  ang180 = mod(ang, 360.0) - (sign(1.0,ang)-1.0)*180.0
  ang180 = ang180 - 180.0*(sign(1.0, ang180-180.0)+1.0)

end function ang180

! === FUNCTION ang360 ======================================================== !
! Maps arbitrary angle [deg] to [0, 360) [deg].
!
! 2001-08-15, Oyvind.Breivik@dnmi.no
! 2021-06-08, Veronica.Berglyd.Olsen@met.no : Moved from sphere.f
! ============================================================================ !
real pure function ang360(ang)

  implicit none
  real, intent(in) :: ang ! [deg]

  ang360 = mod(ang, 360.0) - (sign(1.0,ang)-1.0)*180.0

end function ang360

end module sphere_vec
