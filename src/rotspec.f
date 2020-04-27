      PROGRAM rotspec
!
! Turn spectra from one rotated spherical projection to another.
! Note: Handles only rotated spherical projections.
!
! USAGE: rotspec inspec.seq outspec.seq xcen ycen [xcen2 ycen2]
!
! Input  file   - inspec.seq:  sequential binary 2D spectra - DNMI format
! Output file   - outspec.seq: same format as input file, spectra turned to 
!                  desired projection
! xcen ycen     - lon and lat of center pos in old rotated spherical grid,
! [xcen2 ycen2] - lon and lat of center pos in new rotated spherical grid,
!                  default to xcen2=0.0, ycen2=0.0 (unrotated grid) if omitted
!
! Requires:
!    sphere.f
!    rotspher.f
!    specio.f
!
! Compile IBM:
! f90 -qfixed -O2 -C rotspec.f sphere.f rotspher.f specio.f -o rotspec
! Compile Linux:
! gfortran -O2 -C rotspec.f sphere.f rotspher.f specio.f -o rotspec
!
! 2006m01d23, Oyvind.Breivik@met.no
! 2007m10d15, changed from wamspecio to specio.
! 2009-05-27, added support for multiple analysis times
! 2009-05-28, bug fix: handles single time

! Test:
! % rotspec spk00.dat spk00rot.dat 0.0 0.0

      IMPLICIT NONE

      ! Functions called
      REAL degrees
      REAL ang360
      INTEGER iargc

      ! Constants
      INTEGER MAXANG, MAXFRE, MAXARR
      INTEGER IU25, IU26
      REAL XCEN0, YCEN0

      PARAMETER ( MAXANG=50, MAXFRE=51, MAXARR=4000 )
      PARAMETER ( IU25=25, IU26=26 )
      PARAMETER ( XCEN0=0.0, YCEN0=0.0 )

      ! Vars
      INTEGER fchh, fchh0
      INTEGER immdd0, immdd, ihhmi0, ihhmi
      INTEGER ntime, nang, nfre, nspec, narg
      INTEGER dk, i, k, k2, kplus, l, m, n
      INTEGER ierr
      REAL pi, rad, deg
      REAL lonrad, latrad, rlon, rlat
      REAL xcen, ycen, xcen2, ycen2
      REAL thw, thq
      REAL delth, angle, frac, delta

      ! Arrays
      INTEGER ident(20), identnew(20), gpdummy(6), itime(5)
      REAL slon(MAXARR), slat(MAXARR)
      REAL delang(MAXARR), delang2(MAXARR)
      REAL sp(MAXANG,MAXFRE), spr(MAXANG,MAXFRE)

      ! Strings
      CHARACTER*256 fin, fout, str

      ! Initialize
      xcen=XCEN0
      ycen=YCEN0
      ierr=0
      ntime=0
      nspec=0
      pi=4.0*atan(1.0)
      rad=pi/180.0
      deg=180.0/pi

      ! Command line
      narg=iargc()
      if (narg < 3) then
         write (*,*) 
     2   "USAGE: rotspec inspec.seq outspec.seq xcen ycen [xcen2 ycen2]"
         write (*,*)
     2   " Rotate 2D spectra between two rotated spherical grids"
         write (*,*) 
     2   " inspec.seq: binary sequential input file with 2D spectra"
         write (*,*) 
     2   " outspec.seq: output file w/spectra rotated to new projection"
         write (*,*) 
     2   " xcen ycen: Lon and lat of center point [deg] of old grid,"
         write (*,*) 
     2   " xcen2 ycen2: Center point of new rot. spherical grid,"
         write (*,*) 
     2   "  xcen2=0.0, ycen2=0.0 [default] converts to std lon/lat"
         stop
      endif

      ! Open binary sequential file with 2D spectra
      call getarg(1,fin)
      open (IU25, file=fin, access="sequential", status="old",
     2      form="unformatted", err=1025)

      ! Output file
      call getarg(2,fout)
      open (IU26, file=fout, access="sequential", status="unknown",
     2      form="unformatted", err=1026)

      ! Read grid parameters xcen, ycen, xcen2, ycen1
       call getarg(3,str)
       read (str,*) xcen 
       call getarg(4,str)
       read (str,*) ycen 
       xcen2 = XCEN0
       ycen2 = YCEN0 ! Defaults to unrotated grid
       if (narg > 4) then
          call getarg(5,str)
          read (str,*) xcen2
          call getarg(6,str)
          read (str,*) ycen2
       endif


!!! Count spectral grid points !!!

      ! Read first spectrum
      call getspec(IU25, ident, MAXANG, MAXFRE, sp, ierr)
      if (ierr /= 0) goto 1010

      nang = ident(10)
      nfre = ident(11)
      delth = 360.0/real(nang)
      fchh0 = ident(4) ! first forecast time
      fchh = fchh0
      ihhmi0 = ident(3)
      ihhmi = ihhmi0
      immdd0 = ident(2)
      immdd = immdd0
      
      do while ( ierr==0 .AND. fchh==fchh0 .AND. ihhmi==ihhmi0 .AND. 
     +           immdd==immdd0 )
         nspec = nspec+1
         slat(nspec) = degrees(ident(5),ident(9))
         slon(nspec) = degrees(ident(6),ident(17))
         latrad = slat(nspec)*rad
         lonrad = slon(nspec)*rad

         ! Compute coordinates in old grid
         call sphrot(lonrad,latrad,rlon,rlat,1,xcen*rad,ycen*rad,1)
         ! Radians to degrees
         rlon=rlon*deg
         rlat=rlat*deg
         ! Local distortion angle in old grid [deg]
         delang(nspec) = real(ident(14))/100.0
         call ancorr(delta,rlon,rlat,xcen,ycen)

         ! Compute coordinates in new grid
         call sphrot(lonrad,latrad,rlon,rlat,1,xcen2*rad,ycen2*rad,1)
         ! Radians to degrees
         rlon=rlon*deg
         rlat=rlat*deg
         ! Local distortion angle in new grid [deg]
         call ancorr(delta,rlon,rlat,xcen2,ycen2)
         delang2(nspec) = delta

         call getspec(IU25, ident, MAXANG, MAXFRE, sp, ierr)
         !if (ierr /= 0) goto 1010
         fchh = ident(4)
         ihhmi = ident(3)
         immdd = ident(2)

      enddo ! while

      ! Rewind spectral input file
      ierr = 0
      rewind IU25


!!! Time loop !!!

      do while ( ierr==0 )

         ntime = ntime+1

         ! Loop over spectral locations
         do l = 1, nspec

            ! Read spectrum
            call getspec(IU25, ident, MAXANG, MAXFRE, sp, ierr)
            if (ierr /= 0) goto 1000

            ! Compute angle between old and new grid at location l
            angle = ang360(delang2(l)-delang(l))

            ! Compute gap measured in directional bins
            dk = int(angle/delth)
            frac = 1.0-mod(angle,delth)/delth

            ! Rotate to new grid
            do k2 = 1, nang
               !k = modulo(k2+dk-1, nang) + 1
               k = mod(k2+dk-1, nang) + 1
               !kplus = modulo(k, nang) + 1
               kplus = mod(k, nang) + 1
               do m = 1, nfre
                  spr(k2,m) = frac*sp(k,m)+(1.0-frac)*sp(kplus,m)
               enddo
            enddo

            ! Prepare new header
            do i = 1, 20
               identnew(i) = ident(i)
            enddo

            ! Wind direction rel to new rotated N (going to)
            thw = real(ident(8))/10.0
            thw = ang360(thw-angle)
            identnew(8) = nint(thw*10.0) ! [deg/10]

            ! Angle correction from rotated meridian to true meridian
            identnew(14) = nint(delang2(l)*100.0)

            ! Wave direction rel to new rotated N (going to)
            thq = real(ident(18))/10.0
            thq = ang360(thq-angle)
            identnew(18) = nint(thq*10.0) ! [deg/10]

            ! Automatic scaling
            identnew(20) = -32767

            ! Store spectrum and scalars
            call putspec(IU26, identnew, MAXANG, MAXFRE, spr, ierr)

         enddo ! l

      enddo ! while ntime

!!! End time loop !!!


      ! Summary
1000  continue
      ! Adjust time counter
      ntime = ntime-1
      if (nspec==1) then
         write (*,*) "rotspec: wrote one spectrum at",ntime,
     2               " times"
      else if (nspec>1) then
         write (*,*) "rotspec: wrote", nspec, " spectra at",ntime,
     2               " times"
      endif
         
      goto 1030

      ! Error handling
1010  continue
      write (*,*) "ERROR: rotspec: Error in getspec, ierr == ", ierr
      goto 1030
1025  continue
      write (*,*) "ERROR: rotspec: Open error unit ", IU25
      goto 1030
1026  continue
      write (*,*) "ERROR: rotspec: Open error unit ", IU26
      goto 1030

      ! Close files
1030  continue
      close (IU25)
      close (IU26)

      end
