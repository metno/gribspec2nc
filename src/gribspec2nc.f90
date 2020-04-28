! ----------------------------------------------------------------------
!
PROGRAM gribspec2nc
!
!     VERSION:
!     ---------

!     Oyvind Breivik, Norwegian Meteorological Institute, 2015-02-25
!     oyvind.breivik@met.no (OeB)
!     2016-08-15: Added option -s to issue error message or warning when encountering irregularities in the time stepping
!     2017-06-08: Fixed bug in computation of time step. OeB.
!
!
!     PURPOSE
!     -------
!     Decodes GRIB spectra stored as GRIB parameter 251.
!     The results will be saved as a NetCDF file containing
!     spectra in selected locations provided as input.
!     The spectra are interpolated to these locations if possible.
!
!      USAGE: gribspec2nc [-x] [-d timestep_hh][-i infile] [-l speclist.inp] [-o outfile] [-t itest] [-w weightfile]
!               -x: compute integrated parameters
!               -d timestep_hh: specify temporal resolution, default is 1 h
!               -i infile: GRIB input, default is input_spectra
!               -l speclist.inp: list of spectral locations to interpolate to,
!                           default is speclist.inp
!               -o outfile: NetCDF output, default is output_spectra.nc
!               -w weightfile: Use precomputed interpolation weights, default is no
!               -s: issue a warning but generate NetCDF file even if there are missing GRIB fields (default is to issue error and abort)
!               -t: debug reporting:
!                   itest: 0 (no diagnostics)
!                   itest: 1 (some diagnostics)
!                   itest: 2 (some more diagnostics)
!                   itest: 3 (all diagnostics)
!
!     INPUT FILE REQUIREMENT :
!     -----------------------
!     The input file can only contain one parameter (251) obtained
!     for all directions and frequencies in the default order. Namely, the
!     mars request should be done with DIRECTION=1/TO/nang, and
!     FREQUENCY=1/TO/nfre. At this time nfre=30 and nang=24 for global model
!     and for mediterranean data.
!
!     LIBRARY : FORTRAN90, GRIBAPI
!     -------
!     COMPILE:
!     export NETCDF4_DIR="/opt/netcdf-fortran-4.2"
!     export NETCDF_INCLUDE="-I$NETCDF4_DIR/include"
!     export NETCDF_LIB="-L$NETCDF4_DIR/lib -lnetcdff"
!     gfortran gribspec2nc.f90 sphere.f netcdf_metno_spec.f90 $NETCDF_INCLUDE $NETCDF_LIB $GRIB_API_INCLUDE $GRIB_API_LIB -o gribspec2nc
!
!     export GRIB_API_LIB="-L/opt/grib_api-1.13.0/lib -Wl,-rpath,/opt/grib_api-1.13.0/ -lgrib_api_f90 -lgrib_api -lm"
!     export NETCDF4_DIR="/usr/lib/x86_64-linux-gnu/"
!     export NETCDF_INCLUDE="-I/usr/include"
!     export NETCDF_LIB="-L$NETCDF4_DIR -lnetcdff"
!     gfortran gribspec2nc.f90 sphere.f netcdf_metno_spec.f90 $NETCDF_INCLUDE $NETCDF_LIB $GRIB_API_INCLUDE $GRIB_API_LIB -o gribspec2nc
!
!     ECMWF:
!     module load netcdf4
!     gfortran gribspec2nc.f90 sphere.f netcdf_metno_spec.f90 $NETCDF_INCLUDE $NETCDF_LIB $GRIB_API_INCLUDE $GRIB_API_LIB -o gribspec2nc
!
!     Run:
!     gribspec2nc -i input_spectra -l speclist.inp -o output_spectra.nc
!     gribspec2nc -x -i ../Src/ifs_spectra_an.20150802.grb -l ifs_spectra_list.inp -w ifs_weights.bin -o ifs_spectra_an.20150802.nc
!     Analysis file (covers WAM4 domain):
!     gribspec2nc -d 1 -x -i ~/Data/Mdata/Spec/ifs_spectra_wam4.20160111_an.grb -l invdist_MYWAVE4.inp -w ifs_weights_wam4.bin -o ifs_spectra_an_wam4.nc
!     Forecast file:
!     gribspec2nc -d 1 -x -i ~/Data/Mdata/Spec/ifs_spectra.20160110_fc.grb -l ifs_spectra_list.inp -w ifs_weights.bin -o ifs_spectra_fc.nc
!     ERA-Interim file:
!     gribspec2nc -x -d 6 -i /metno/hindcast3/oyvindb/Sea_of_Japan//erai_soj_spectra_20140101-20140115.grb -l newinvsJapan.inp -o erai_spectra_soj.20140101-20140115.nc
!     gribspec2nc -d 6 -x -i ~/Data/Mdata/Spec/erai_soj_spectra_201401.grb -l newinvsJapan.inp -o ~/Data/Mdata/Spec/erai_spectra_soj.201401.nc
!     COMPILE:
!     gfortran gribspec2nc.f90 sphere.f netcdf_metno_spec.f90 $NETCDF_INCLUDE -L/opt/netcdf-fortran-4.2/lib -lnetcdff $GRIB_API_INCLUDE $GRIB_API_LIB -o gribspec2nc
!
!
!     Run:
!     gribspec2nc -i input_spectra -l speclist.inp -o output_spectra.nc
!     gribspec2nc -x -i ../Src/ifs_spectra_an.20150802.grb -l ifs_spectra_list.inp -w ifs_weights.bin -o ifs_spectra_an.20150802.nc
!     Analysis file (covers WAM4 domain):
!     gribspec2nc -d 1 -x -i ~/Data/Mdata/Spec/ifs_spectra_wam4.20160111_an.grb -l invdist_MYWAVE4.inp -w ifs_weights_wam4.bin -o ifs_spectra_an_wam4.nc
!     Forecast file:
!     gribspec2nc -d 1 -x -i ~/Data/Mdata/Spec/ifs_spectra.20160110_fc.grb -l ifs_spectra_list.inp -w ifs_weights.bin -o ifs_spectra_fc.nc
!
!
!     MARS REQUEST EXAMPLE:
!     -----------------------

! ----------------------------------------------------------------------

   USE GRIB_API
   USE netcdf
   USE netcdf_metno_spec
   use sphere_vec, only : spheredist

   IMPLICIT NONE

   ! Functions called
  !  REAL :: spheredist
   REAL :: ang180
   REAL :: ang360
   INTEGER :: getclo
   INTEGER :: julday

   ! Constants
   REAL, PARAMETER :: ZMISS=-999.0
   REAL, PARAMETER :: PI=3.1415927
   REAL, PARAMETER :: ZPI=2.0*PI
   REAL, PARAMETER :: G = 9.806
   REAL, PARAMETER :: RAD = PI/180.0
   REAL, PARAMETER :: DEG = 180.0/PI
   REAL, PARAMETER :: TOL=0.1
   REAL, PARAMETER :: POW=2.0
   REAL, PARAMETER :: DISTMAX=40000.0

   INTEGER, PARAMETER :: NOUT_S = 6        !! NUMBER OF OUTPUT PARAMS.
   INTEGER, PARAMETER :: IFSPECLIST = 21
   INTEGER, PARAMETER :: IFWEIGHTS = 23

   ! Vars
   INTEGER :: morarg, iu06
   INTEGER :: i1, itest, lfile
   INTEGER :: itabpar, iparam, itable, irgg, iscan
   INTEGER :: iyyyymmdd, ihhmm, ihh, imm, idd, iyyyy, imi, ihha, ihh0
   INTEGER :: jday0, jdaya
   INTEGER :: ngx, ngy, l, k, i, j, m, jsn, istart, ij, kmax, mmax
   INTEGER :: nang, nfre, nfrang, kk, mm
   INTEGER :: idirscaling, ifrescaling
   INTEGER :: ilon, ilat, ierr
   INTEGER :: nwish, neighbours, miss
   INTEGER :: getcla, ioptval
   INTEGER :: ifile, iret, igrib, ios
   INTEGER :: numberofvalues, nbval
   INTEGER :: ngxloc
   INTEGER :: plpresent, nb_pl
   INTEGER :: istep=0, itstep, ifchh
   INTEGER :: ideldo=3600 ! Default

   REAL :: one, ztheta, zfre
   REAL :: zdello
   REAL :: yfrst, ylast, amonop, amosop, amowep, amoeap
   REAL :: xdello, xdella
   REAL :: dello
   REAL :: step
   REAL :: dtresh, fchh
   REAL :: delth, delt25, th0, denom, etot, emax, temp, xx
   REAL :: co, co1
   CHARACTER :: c

   LOGICAL :: llexist, llfrst, lfirsterror=.TRUE.
   LOGICAL :: llscanns
   LOGICAL :: lperiodic
   LOGICAL :: llresetmissing
   LOGICAL :: lfirst = .TRUE.
   LOGICAL :: lstrict = .TRUE.
   LOGICAL :: lnullspec = .FALSE.
   LOGICAL :: lsaveweights = .FALSE.
   LOGICAL :: lextras = .FALSE.
   LOGICAL :: llexistweights = .FALSE.

   ! Arrays

   LOGICAL, DIMENSION(NOUT_S) :: CFLAG_S    !! COMPUTATION FLAG.

   CHARACTER*1  :: CLOPTLET
   CHARACTER*3  :: CLL1
   CHARACTER*16 :: CLOPTS
   CHARACTER*8  :: CSTEPUNITS
   CHARACTER*8  :: CSTEPTYPE
   CHARACTER*12 :: CGRIDTYPE
   CHARACTER*12 :: CLFMT
   CHARACTER*128:: clarg, fnamein, fnameout, fspeclist, fweights
   CHARACTER*256:: str
   CHARACTER (LEN=14) :: cdatea, cdate0="first", cdatef    !! START DATE OF PRINT OUPUT  (YYYYMMDDHHMMSS).

   ! Allocatables
   INTEGER, DIMENSION(:), ALLOCATABLE :: pl
   INTEGER, DIMENSION(:), ALLOCATABLE :: kdomrgg
   INTEGER, ALLOCATABLE :: idx(:,:)

   REAL(KIND=8), ALLOCATABLE :: values(:)
   REAL, ALLOCATABLE :: spec(:)
   REAL, ALLOCATABLE :: scfr(:), scth(:), esumth(:), esumfr(:)
   REAL, ALLOCATABLE :: dfim(:), fr(:), theta(:)
   REAL, ALLOCATABLE :: xlon(:), xlat(:), lon(:), lat(:)
   REAL, ALLOCATABLE :: dist(:), sumw(:)
   REAL, ALLOCATABLE :: w(:,:), distmin(:,:)
   REAL, ALLOCATABLE :: spw(:,:)
   REAL, ALLOCATABLE :: hm0(:,:), pdir(:,:), tpeak(:,:)


   ! Init
   DATA CLOPTS/'i;d;l;o;s;t;w;x;'/
   DATA CFLAG_S/.TRUE.,.FALSE.,.FALSE.,.TRUE.,.TRUE.,.TRUE./

! ----------------------------------------------------------------------

!*    INITIAL VALUES SET AND CRACK COMMAND LINE.
!     -----------------------------------------

   fnamein='input_spectra'
   fnameout='output_spectra.nc'
   fspeclist='speclist.inp'
   fweights=' '
   itest=0
   iu06=6

   CMDLINE: DO
      IOPTVAL=GETCLO(CLOPTS,CLARG)
      IF (IOPTVAL .LE. 0 )  THEN
         EXIT CMDLINE
      ENDIF
      CLOPTLET=CHAR(IOPTVAL)
      ! GETS VARIABLE ARGUMENT FOR OPTION
      morarg=getcla(clarg)
      IF (morarg /= 0) THEN
         IF ( cloptlet == 'i' ) THEN
            fnamein=clarg
         ELSEIF ( cloptlet == 'd' ) THEN
            READ (clarg, *) ideldo
            ideldo = ideldo*3600
         ELSEIF ( cloptlet == 'l' ) THEN
            fspeclist=clarg
         ELSEIF ( cloptlet == 'o' ) THEN
            fnameout=clarg
         ELSEIF ( cloptlet == 't' ) THEN
            i1 = len_trim(clarg)
            WRITE (CLL1,'(I3)') I1
            clfmt = '(I'//CLL1//')'
            READ (CLARG(1:I1),FMT=CLFMT) ITEST
         ELSEIF ( cloptlet == 'w' ) THEN
            fweights=clarg
            lsaveweights = .TRUE.
         ENDIF
      ELSE
         IF ( cloptlet == 'x' ) THEN
            lextras = .TRUE.
         ELSEIF ( cloptlet == 's' ) THEN
            lstrict = .FALSE.
         ENDIF
      ENDIF

   ENDDO CMDLINE


!*    OPEN ASCII SPECTRAL LOCATION FILE AND READ WISH LIST
!     ----------------------------------------------------
   LFILE=0
   IF (fspeclist.NE. ' ') LFILE=LEN_TRIM(fspeclist)
   OPEN (IFSPECLIST, FILE=fspeclist(1:LFILE), form="formatted")

   ! Count number of desired locations
   READ (IFSPECLIST, "(a)", iostat=ios) str
   nwish = 0
   DO WHILE ( .NOT. is_iostat_end(ios) )
      READ (IFSPECLIST, "(a)", iostat=ios) str
      ! List may contain comments, first character must be c, *, ! or #
      IF ( c/='c' .AND. c/='*' .AND. c/='!' .AND. c/='#') THEN
         nwish = nwish+1
      ENDIF
   ENDDO
   ! Ignore first line which contains metadata
   nwish = nwish-1
   ALLOCATE (lon(nwish))
   ALLOCATE (lat(nwish))

   ! Loop until end of wish list
   REWIND (IFSPECLIST)
   READ (IFSPECLIST, "(a)", iostat=ios) str
   j = 0
   lnullspec = .FALSE.
   DO WHILE ( .NOT. is_iostat_end(ios) )
      c = str(1:1)
      ! List may contain comments, first character must be c, *, ! or #
      IF ( c/='c' .AND. c/='*' .AND. c/='!' .AND. c/='#') THEN
         ! First line contains
         ! dtresh - treshold [km] distance to neighbours and
         ! neighbours - max number of neighbours used for interpolation
         IF (j == 0) THEN
            READ (str,*) dtresh, neighbours
            IF (dtresh < 0.0) lnullspec = .TRUE.
            dtresh = abs(dtresh)
         ELSE
            ! Following lines contain lat, lon [deg] of desired locations
            READ (str,*) lat(j), lon(j)
            lon(j) = ang180(lon(j))
         ENDIF
         j = j+1
      ENDIF
      READ (IFSPECLIST, "(a)", iostat=ios) str
   ENDDO
   CLOSE (IFSPECLIST)

   ALLOCATE (sumw(nwish))
   ALLOCATE (idx(neighbours,nwish))
   ALLOCATE (w(neighbours,nwish))
   ALLOCATE (distmin(neighbours,nwish))
   ALLOCATE (hm0(nwish,1))
   ALLOCATE (pdir(nwish,1))
   ALLOCATE (tpeak(nwish,1))

!*    OPEN DATA FILE
!     --------------
   lfile=0
   llexist=.FALSE.
   IF (fnamein /= ' ') lfile=LEN_TRIM(fnamein)
   INQUIRE (FILE=fnamein(1:lfile), EXIST=llexist)
   IF (llexist) THEN
      CALL grib_open_file(ifile,fnamein(1:lfile),'r')
   ELSE
      WRITE(*,*)'****************************'
      WRITE(*,*)'*                          *'
      WRITE(*,*)'*GRIB DATA NOT FOUND IN *'
      WRITE(*,*)  FNAMEIN
      WRITE(*,*)'*PROGRAM WILL ABORT        *'
      WRITE(*,*)'*                          *'
      WRITE(*,*)'****************************'
      CALL ABORT
   ENDIF

!  GET FIRST DATA FILE
!  LOOP ON ALL MESSAGES IN INPUT FILE

   igrib=-99
   CALL grib_new_from_file(ifile,igrib,iret)

   LOOP: DO WHILE (iret /= grib_end_of_file)

      !* DETERMINE DATA FIELD CHARACTERISTICS

      CALL grib_get(igrib, 'paramId', itabpar)
      itable=itabpar/1000
      iparam=itabpar-itable*1000
      IF (itest>0) WRITE (*,*) ' THE INPUT PARAMETER IS ', iparam

      CALL grib_get(igrib,'gridType', cgridtype)
      IF (cgridtype(1:7) == 'regular') THEN
         irgg=0
      ELSEIF (cgridtype(1:7) == 'reduced') THEN
         irgg=1
      ELSE
         WRITE(IU06,*) '***********************************'
         WRITE(IU06,*) '*  GRID TYPE NOT RECOGNIZED !!!'
         WRITE(IU06,*) '   gridType = ', CGRIDTYPE
         WRITE(IU06,*) '***********************************'
         CALL ABORT
      ENDIF

      CALL grib_get(igrib, 'iScansNegatively', iscan)
      IF (iscan==0) THEN
         llscanns=.TRUE.
      ELSEIF (iscan==64 .OR. iscan==2) THEN
         llscanns=.FALSE.
      ELSE
         WRITE (IU06,*) '***********************************'
         WRITE (IU06,*) '*  SCANNING MODE NOT RECOGNIZED !!!'
         WRITE (IU06,*) ' ISCAN = ', ISCAN
         WRITE (IU06,*) '***********************************'
         CALL ABORT
      ENDIF

      CALL grib_get(igrib,'latitudeOfFirstGridPointInDegrees',yfrst)
      CALL grib_get(igrib,'latitudeOfLastGridPointInDegrees',ylast)

      IF (llscanns) THEN
         amonop = yfrst
         amosop = ylast
      ELSE
         amonop = ylast
         amosop = yfrst
      ENDIF

      CALL grib_get(igrib,'longitudeOfFirstGridPointInDegrees',amowep)
      CALL grib_get(igrib,'longitudeOfLastGridPointInDegrees',amoeap)

      CALL grib_get(igrib,'jDirectionIncrementInDegrees',xdella)
      CALL grib_get(igrib,'iDirectionIncrementInDegrees',xdello)

      CALL grib_get(igrib,'Ny',ngy)

      IF (lfirst) ALLOCATE (kdomrgg(ngy))
      kdomrgg=0
      IF (irgg==0) THEN
         CALL grib_get(igrib,'Nx',ngx)
         kdomrgg = ngx
      ELSE
         CALL grib_get(igrib,'PLPresent',plpresent)
         IF (plpresent == 1) THEN
            CALL grib_get_size(igrib,'pl',nb_pl)
            ALLOCATE (pl(nb_pl))
            CALL grib_get(igrib,'pl',pl)
         ELSE
            WRITE (*,*) 'NUMBER OF POINTS PER LATITUDE MISSING !!!'
            CALL abort
         ENDIF

         istart=1
         DO WHILE(pl(istart) == 0 .AND. istart < nb_pl)
            istart = istart+1
         ENDDO
         istart=istart-1
         ngx = 0
         DO j=1,ngy-istart
            kdomrgg(j) = pl(j+istart)
            ngx = MAX(ngx,kdomrgg(j))
         ENDDO
         DEALLOCATE (pl)
      ENDIF

      CALL adjust(amowep, amoeap)
      lperiodic = .false.
      dello = (amoeap-amowep)/MAX(1,ngx-1)
      IF (amoeap-amowep+1.5*dello>=360.0) lperiodic=.TRUE.

      CALL grib_get(igrib, 'dataDate', iyyyymmdd)
      CALL grib_get(igrib, 'time', ihhmm)
      ihh = ihhmm/100
      imi = mod(ihhmm,100)
      iyyyy = iyyyymmdd/10000
      imm = mod(iyyyymmdd,10000)
      idd = mod(imm,100)
      imm = imm/100

      CALL grib_get(igrib, 'stepType', csteptype)

 !       FORECAST STEP (defined here in hours)
      IF (csteptype(1:7) == 'instant') THEN
         CALL grib_set(igrib,'stepUnits','h')
         CALL grib_get(igrib,'step',step)
         CALL grib_get(igrib,'stepRange',fchh)
      ELSE
         WRITE(*,*) 'UNKNOWN DEFINITION OF FORECAST STEP TYPE !!!'
         WRITE(*,*) 'stepType = ',CSTEPTYPE
         CALL ABORT
      ENDIF
      ! hours
      ifchh = nint(fchh)

      WRITE (cdatea, "(i8,i4.4,i2.2)") iyyyymmdd, ihhmm, 0
      call incdate(cdatea,ifchh*3600)
      jdaya = julday(iyyyy,imm,idd)+(ihh+ifchh)/24
      ihha = MOD(ihh+ifchh,24)

      ! Keep initial date (not necessarily analysis time since fchh may be nonzero)
      IF (cdate0 == "first") THEN
         WRITE (cdate0, "(i8,i4.4,i2.2)") iyyyymmdd, ihhmm, 0
         CALL incdate(cdate0,ifchh*3600)
         cdatef = cdate0
         ! Keep Julian day and hour of first record in GRIB file
         jday0 = jdaya
         ihh0 = ihha
      ENDIF
      print *, "jdaya jday0 ihha ihh0", jdaya,jday0,ihha,ihh0
      itstep = NINT(3600*REAL(24*(jdaya-jday0)+ihha-ihh0)/ideldo)


      ! Compare cdatea, date inferred from GRIB file to cdatef, the date inferred assuming a constant time step ideldo
      IF (cdatea /= cdatef) THEN
          IF (lstrict) THEN
             print *, "Gribspec2nc: Error: irregular timestepping or wrong user-specified time step at cdatea = ", cdatea
             print *, "Gribspec2nc: Aborting"
             STOP
          ELSE
             print *, "Gribspec2nc: Warning: irregular timestepping or wrong user-specified time step at cdatea = ", cdatea
          ENDIF
          print *, "istep, itstep ", istep, itstep
      ENDIF
      CALL incdate(cdatef, ideldo)

      IF (iparam==250) THEN
         ! obsolete parameter 250 not handled
         write(*,*) 'not ready yet for 250 !!!!'
         call abort
      ELSEIF (iparam==251) THEN
         CALL grib_get(igrib,'numberOfDirections',nang)
         CALL grib_get(igrib,'numberOfFrequencies',nfre)
      ELSE
         WRITE(*,*) 'THE INPUT GRIB PARAMETER IS NOT 250 OR 251 BUT', IPARAM
         WRITE(*,*) 'WHICH IS NOT A WAVE SPECTRUM PARAMETER !!!'
         WRITE(*,*) 'PROGRAM WILL ABORT'
         CALL ABORT
      ENDIF

      nfrang=nang*nfre

      IF (lfirst) THEN
         ALLOCATE (fr(nfre))
         ALLOCATE (scfr(nfre))
         ALLOCATE (esumth(nfre))
         ALLOCATE (theta(nang))
         ALLOCATE (esumfr(0:nang+1))
         ALLOCATE (scth(nang))
         ALLOCATE (spw(nang,nfre))
      ENDIF

      delth = ZPI/nang

 !       DECODE INPUT GRIB DATA
 !       ----------------------

      IF (iparam == 250) THEN

         WRITE (*,*) ' not yet !!! '

      ELSEIF (iparam == 251) THEN

         ! GET THE SIZE OF THE VALUES ARRAY
         CALL grib_get_size(igrib,'values',numberofvalues)

         IF (lfirst) THEN
            ALLOCATE (values(numberofvalues))
            ALLOCATE (dist(numberofvalues))
            ALLOCATE (spec(numberofvalues*nfrang))
         ENDIF

         DO m = 1,nfre
            DO k = 1,nang
   !           CHECK ON NUMBER OF VALUES
               CALL GRIB_GET_SIZE(IGRIB,'values',NBVAL)
               IF (nbval /= numberofvalues) THEN
                  WRITE (*,*) '*************************************'
                  WRITE (*,*) 'THE NUMBER OF VALUES SHOULD ALWAYS BE'
                  WRITE (*,*) 'NUMBEROFVALUES = ',NUMBEROFVALUES
                  WRITE (*,*) 'IT IS NOW ',NBVAL
                  WRITE (*,*) 'PROBLEM !!!! ABORTING NOW !'
                  WRITE (*,*) '*************************************'
                  !CALL ABORT
                  STOP
               ENDIF

    !          SET THE MISSING DATA INDICATOR
               CALL grib_set(igrib, 'missingValue', zmiss)
               values=zmiss

    !          GET DATA VALUES
               CALL grib_get(igrib,'values',values)

    !          DETERMINE DATA FIELD CHARACTERISTICS

               CALL grib_get(igrib,'directionNumber',kk)
               CALL grib_get(igrib,'directionScalingFactor',idirscaling)
               CALL grib_get(igrib,'scaledDirections',scth)
               theta(kk) = scth(kk)/idirscaling
               ! First bin in degrees
               th0 = theta(1)

               CALL grib_get(igrib,'frequencyNumber',mm)
               CALL grib_get(igrib,'frequencyScalingFactor',ifrescaling)
               CALL grib_get(igrib,'scaledFrequencies',scfr)
               fr(mm) = scfr(mm)/ifrescaling
               delt25 = fr(nfre)/4.0*delth

               IF (kk/=k .AND. mm/=m .AND. lfirsterror) THEN
                  WRITE (*,*) '* Warning: jumbled directions and freqencies.'
                  WRITE (*,*) '* Please check NetCDF spectra.'
                  lfirsterror = .FALSE.
                  !CALL ABORT
               ENDIF

    !          there is currently a bug in the grib_api that sets all values to 0
    !          when they are all missing. Reset them to missing (it will hopefully be
    !          corrected soon).
               llresetmissing=.true.
               DO ij=1,numberofvalues
                  IF (values(ij)/=0.0) THEN
                     llresetmissing=.false.
                     EXIT
                  ENDIF
               ENDDO
               IF (llresetmissing) values=zmiss

               DO ij = 1, numberofvalues
                  IF (VALUES(IJ) /= ZMISS) THEN
                     spec(kk+(mm-1)*nang+(ij-1)*nfrang) = 10.**values(ij)
                  ELSE
                     spec(kk+(mm-1)*nang+(ij-1)*nfrang) = 0.0
                  ENDIF
               ENDDO

    !          GET NEXT FIELD

               IF ( .NOT. (k==nang .AND. m==nfre) ) THEN
                  CALL grib_release(igrib)
                  igrib=-99
                  CALL grib_new_from_file(ifile, igrib, iret)
               ENDIF
            ENDDO ! nang
         ENDDO ! nfre
      ENDIF ! param 250 or 251

     ! Initialize weights etc. Only once

      IF (lfirst) THEN
         !*    OPEN NetCDF OUTPUT FILE
         !     ---------------------------------------
         CALL nc_open_specfile(fnameout, cdatea, cflag_s, lon, lat, fr, theta, ideldo)
         print *, "cdatea = ",cdatea
         WRITE (iu06,*) ' +++ NetCDF file opened ', trim(fnameout)

        ! Delta-frequency array
         ALLOCATE (dfim(nfre))
         ALLOCATE (xlon(numberofvalues))
         ALLOCATE (xlat(numberofvalues))

         co = fr(2)/fr(1)
         co1 = 0.5*(co-1.0)*delth
         dfim(1) = co1*fr(1)
         DO m = 2, nfre-1
           dfim(m) = co1*(fr(m)+fr(m-1))
         ENDDO
         dfim(nfre) = co1*fr(nfre-1)

        ! Find lat and lon of all grid points
         ilon = 0
         ilat = 1
         DO ij = 1, numberofvalues
           ilon = ilon+1
           ! Regular grid ...
           IF (irgg==0) THEN
             IF (ilon>ngx) THEN
               ilon = 1
               ilat = ilat+1
             ENDIF
             zdello = xdello
             ! ... or only wet points
             ELSE
               ngxloc = kdomrgg(ilat)
               IF (ilon > ngxloc) then
                 ilon = 1
                 ilat = ilat+1
                 IF (lperiodic) THEN
                   zdello = 360.0/ngxloc
                 ELSE
                   zdello = (amoeap-amowep)/(ngxloc-1)
                 ENDIF
              ENDIF
           ENDIF

          ! Compute longitude and latitude
           xlon(ij)=ang180(AMOWEP+(ilon-1)*zdello)
           xlat(ij)=amonop-(ilat-1)*xdella
         ENDDO ! ij = 1, numberofvalues

       ! Compute weights based on distances from wishlist to all grid points

         ! Read previously generated weight file?
         IF (lsaveweights) THEN
           IF (fweights /= ' ') lfile=LEN_TRIM(fweights)
           INQUIRE (FILE=fweights(1:lfile), EXIST=llexistweights)
           IF (llexistweights) THEN
              OPEN (IFWEIGHTS, FILE=fweights(1:lfile), form="unformatted", status="old")
              READ (IFWEIGHTS) sumw, idx, w, distmin
              WRITE (*,*) "Read weights from file ", fweights
              CLOSE (IFWEIGHTS)
           ENDIF
         ENDIF

         IF (.NOT. llexistweights) THEN
            ! Loop over wishlist
            miss = 0
            DO j = 1, nwish

              ! Loop over grid points
              call spheredist(xlon(:), xlat(:), lon(j), lat(j), dist(:), numberofvalues)
              dist = dist/1000.0
              do ij= 1,numberofvalues
                if(dist(ij) < TOL) then
                  dist(ij) = TOL
                end if
              end do

               ! Find nearest neighbours
               DO i = 1, neighbours
                  distmin(i,j) = DISTMAX

                  ! Loop over grid points
                  DO ij = 1, numberofvalues
                     IF (dist(ij) < distmin(i,j)) THEN
                        distmin(i,j) = dist(ij)
                        idx(i,j) = ij
                     ENDIF
                  ENDDO ! ij = 1, numberofvalues

                  ! Remove nearest spectral location from next search
                  dist(idx(i,j)) = DISTMAX
               ENDDO ! i = 1, neighbours

               ! Compute weights
               i = 1
               sumw(j) = 0.0001
               DO WHILE ((i<=neighbours) .AND. (distmin(i,j)<=dtresh))
                  w(i,j) = distmin(i,j)**(-POW)
                  sumw(j) = sumw(j) + w(i,j)
                  i = i+1
               ENDDO ! i; end while

               ! Flag loners
               IF (distmin(1,j) > dtresh) THEN
                  miss = miss+1
                  IF (itest>0) THEN
                     WRITE (*,"(a,i5,a,f11.6,a,f12.6)") "WARNING: No spectra within range for pos",j," at lat",lat(j),", lon",lon(j)
                  ENDIF
               ENDIF
            ENDDO ! j = 1, nwish

            IF (lsaveweights) THEN
               OPEN (IFWEIGHTS, FILE=fweights(1:LFILE), form="unformatted", status="unknown")
               WRITE (IFWEIGHTS) sumw, idx, w, distmin
               WRITE (*,*) "Saved weights to file ", fweights
               CLOSE (IFWEIGHTS)
            ENDIF ! lsaveweights
         ENDIF ! .NOT. llexistweights

         lfirst = .FALSE.

      ENDIF ! lfirst

     !!! Interpolate spectra to list of locations

      ! Loop over wish list
      hm0 = 0.0
      pdir = 0.0
      tpeak = 0.0
      DO j = 1, nwish

         ! Interpolate spectrum
         spw(:,:) = 0.0

         ! Loop over nearest neighbours
         i = 1
         DO WHILE ( i <= neighbours .AND. distmin(i,j) <= dtresh )
            ! Grid index
            ij = idx(i,j)

            ! Spectral component weights
            DO m = 1, nfre
               DO k = 1, nang
                  spw(k,m) = spw(k,m) + spec(k+(m-1)*nang+(ij-1)*nfrang)*w(i,j)
               ENDDO
            ENDDO
            i = i+1
         ENDDO ! while i

         ! Divide spectrum by sum of weights and compute scalars of those within range
         IF (distmin(1,j) <= dtresh) THEN

            ! Divide spectrum by sum of weights
            spw(:,:) = spw(:,:)/sumw(j)

            ! Compute integrated parameters?
            IF (lextras) THEN

               ! Compute variance and hm0
               etot=0.0
               ! Loop through all spectral bins
               DO m = 1, nfre
                  temp = 0.0
                  DO k = 1, nang
                     temp = temp+spw(k,m)
                  ENDDO
                  etot = etot+temp*dfim(m)
               ENDDO
               ! Compute tail energy
               etot = etot+delt25*temp
               hm0(j,1) = 4.0*sqrt(etot)

               ! Find peak period
               emax = 0.0
               mmax = 0
               DO m = 1, nfre
                  esumth(m) = 0.0
                  DO k = 1,nang
                     esumth(m) = esumth(m)+spw(k,m)
                  ENDDO
                  IF (esumth(m) > emax) THEN
                     emax = esumth(m)
                     mmax = m
                  ENDIF
               ENDDO

               ! More accurate computation of peak period
               IF (mmax == 0) THEN
                  tpeak(j,1) = 0.0 ! Null spectrum
               ELSEIF (mmax == 1) THEN
                  tpeak(j,1)=1.0/fr(1)
               ELSEIF (mmax == nfre) THEN
                  tpeak(j,1)=1.0/fr(nfre)
                  IF (hm0(j,1)<tpeak(j,1)**2/25.0) tpeak(j,1) = 5.0*sqrt(hm0(j,1))
               ELSE
                  xx = (esumth(mmax+1)-esumth(mmax-1))/ &
                       (esumth(mmax+1)-esumth(mmax)-esumth(mmax)+esumth(mmax-1))
                  tpeak(j,1) = (1.0+xx*(.066+2.04315e-3*(1.0-xx)))/fr(mmax)
               ENDIF

               ! Compute peak direction
               emax = 0.0
               esumfr = 0.0
               kmax = 1
               DO k = 1, nang
                  esumfr(k) = 0.0
                  DO m = 1, nfre
                     esumfr(k) = esumfr(k) + spw(k,m)
                  ENDDO
                  IF (esumfr(k) > emax) THEN
                     emax = esumfr(k)
                     kmax = k
                  ENDIF
               ENDDO

               ! More accurate computation of peak direction
               IF (esumfr(kmax) > 0.0) THEN
                  esumfr(0) = esumfr(nang)
                  esumfr(nang+1) = esumfr(1)
                  denom = esumfr(kmax+1)-esumfr(kmax)*2.0 + esumfr(kmax-1)
                  pdir(j,1) = theta(kmax) - 0.5*DEG*delth*(esumfr(kmax+1) - esumfr(kmax-1))/(denom+0.001)
                  pdir(j,1) = ang360(pdir(j,1) + th0)
                  IF (pdir(j,1) < 0.0) pdir(j,1) = pdir(j,1) + 360.0
                  IF (pdir(j,1) >= 360.0) pdir(j,1) = pdir(j,1) - 360.0
               ENDIF

            ENDIF ! lextras

         ENDIF ! distmin <= dtresh

         ! Write spectrum if location is within range or lnullspec is true
         IF ( (distmin(1,j) <= dtresh) .OR. lnullspec) THEN
            CALL nc_write_spec(1,spw,itstep,j,ideldo)
         ENDIF ! ( (distmin(j,1) <= dtresh) .OR. lnullspec)
      ENDDO ! j = 1, nwish
      IF (lextras) THEN
         CALL nc_write_intpar(4,hm0,itstep,ideldo)
         CALL nc_write_intpar(5,tpeak,itstep,ideldo)
         CALL nc_write_intpar(6,pdir,itstep,ideldo)
      ENDIF ! lextras

      CALL grib_release(igrib)
      igrib=-99

      istep = istep+1
      CALL grib_new_from_file(ifile, igrib, iret)

   ENDDO LOOP

   ! Close GRIB input file
   CALL grib_close_file(ifile)

   ! Close NetCDF output file
   CALL nc_close_specfile

   IF (ALLOCATED(kdomrgg)) DEALLOCATE(kdomrgg)
   IF (ALLOCATED(fr)) DEALLOCATE(fr)
   IF (ALLOCATED(scfr)) DEALLOCATE(scfr)
   IF (ALLOCATED(scth)) DEALLOCATE(scth)
   IF (ALLOCATED(theta)) DEALLOCATE(theta)
   IF (ALLOCATED(values)) DEALLOCATE(values)
   IF (ALLOCATED(spec)) DEALLOCATE(spec)
   IF (ALLOCATED(spw)) DEALLOCATE(spw)
   IF (ALLOCATED(esumth)) DEALLOCATE(esumth)
   IF (ALLOCATED(esumfr)) DEALLOCATE(esumfr)
   IF (ALLOCATED(dist)) DEALLOCATE(dist)

   IF (ALLOCATED(dfim)) DEALLOCATE(dfim)
   IF (ALLOCATED(lon)) DEALLOCATE(lon)
   IF (ALLOCATED(lat)) DEALLOCATE(lat)
   IF (ALLOCATED(xlon)) DEALLOCATE(xlon)
   IF (ALLOCATED(xlat)) DEALLOCATE(xlat)
   IF (ALLOCATED(idx)) DEALLOCATE (idx)
   IF (ALLOCATED(w)) DEALLOCATE (w)
   IF (ALLOCATED(sumw)) DEALLOCATE (sumw)
   IF (ALLOCATED(distmin)) DEALLOCATE (distmin)
   IF (ALLOCATED(hm0)) DEALLOCATE (hm0)
   IF (ALLOCATED(pdir)) DEALLOCATE (pdir)
   IF (ALLOCATED(tpeak)) DEALLOCATE (tpeak)

END PROGRAM gribspec2nc

!#######################################################################
FUNCTION getclo(yaoptions, yaargument)
   INTEGER getclo, getcla
   CHARACTER*   1 yolastarg
   CHARACTER* (*) yaoptions, yaargument
   CHARACTER* 120 arg

   INTEGER here, imorearg, ivarg
   DATA here, imorearg, ivarg, arg / 1, 0, 0, "  " /
   DATA yolastarg / " " /

   arg=' '
   CALL getarg(here,arg)
   iol = len_trim(arg)
   IF (iol == 2 .AND. arg(1:1) == '-' .AND. ivarg == 0 ) THEN
      iol = len_trim(yaoptions)
      DO jl = 1, iol
         getclo = 0
         IF ( yaoptions(jl:jl) .EQ. arg(2:2) ) THEN
            getclo = ichar(arg(2:2))
            IF (yaoptions(jl+1:jl+1) .EQ. ':' ) THEN
               yolastarg=yaoptions(jl:jl)
               ivarg=1
            ENDIF
            EXIT
         ENDIF
      ENDDO
   ELSEIF ( ivarg == 1 ) THEN
      WRITE (*,*) ' option -', yolastarg, ' requires arguments'
      getclo=-1
   ELSEIF (iol == 0) THEN
      getclo=0
   ELSE
      WRITE(*,*) 'illegal option: ',arg(1:iol)
      getclo=-1
   ENDIF
   here = here + 1
   RETURN

   ENTRY getcla(yaargument)
!-->  PRINT*,'-------------getcla--------------'
!-->  PRINT*, 'HERE ins getcla :', here,' options: ', yaoptions

      getcla = 1
      CALL getarg(here,arg)
!-->  PRINT*,' arg in getcla ', arg
      IF ( arg (1:1) .NE. '-' ) THEN
         here = here + 1
         yaargument=arg
      ELSE
         IF (ivarg.EQ.1) THEN
            WRITE(*,*)' refused to take ', arg (1:2) ,' as argument for', ' the option -',yolastarg
            getcla = -1
         ELSE
            getcla = 0
         ENDIF
      ENDIF
      ivarg=0
   !-->  PRINT*,' getcla in getcla ', getcla

   RETURN
END

SUBROUTINE ADJUST(WEST, EAST)

! ----------------------------------------------------------------------

!**** *ADJUST* - ROUTINE TO CORRECT BORDERS OF INTERVALS.

!     H.GUNTHER            ECMWF       04/04/1990

!*    PURPOSE.
!     -------

!       ADJUSTS INTERVAL BORDERS GIVEN IN DEGREE.

!**   INTERFACE.
!     ----------

!       *CALL* *ADJUST (WEST, EAST)*
!          *WEST*    - LEFT INTERVAL BORDER IN DEGREE.
!          *EAST*    - RIGHT INTERVAL BORDER IN DEGREE.

!     METHOD.
!     -------

!       THE INTERVAL BORDERS ARE CHANGED TO FULLFILL:
!         0. .LE. EAST  .AND. EAST .LT. 360. .AND. WEST .LE. EAST

!     EXTERNALS.
!     ----------

!       NONE.

!     REFERENCE.
!     ----------

!       NONE.

! ----------------------------------------------------------------------

!* 1. CORRECT BORDERS.
!     ----------------

   WEST = MOD(WEST+720.,360.)
   EAST = MOD(EAST+720.,360.)
   IF (WEST.GT.EAST) WEST = WEST-360.

   RETURN
END SUBROUTINE ADJUST

! ======================================================================

      SUBROUTINE INCDATE (CDATE, ISHIFT)

! ======================================================================

!**** *INCDATE* - TO UPDATE DATE TIME GROUP

!     J. BIDLOT   FEB 2007    RE-INRODUCING THE OLD WAY WITHOUT
!                             THE NEED FOR ECLIB. ADDING SECONDS

!**   PURPOSE.
!     --------
!       UPDATING DATE TIME GROUP.

!**   INTERFACE.
!     ----------

!       *CALL* *INCDATE (CDATE,ISHIFT)*
!         *CDATE*  CHAR*12 - DATE TIME GROUP (YYYYMMDDHHMM) OR
!                  CHAR*14 - DATE TIME GROUP (YYYYMMDDHHMMSS)
!         *ISHIFT* INTEGER - TIME INCREMENT IN SECONDS, RESTRICTIONS:
!                            SINCE ISHIFT IS AN INTEGER
!                            -MAXINT < ISHIFT <  +MAXINT
!                             WHERE MAXINT=2**(IREP-1)
!                             WHERE IREPR IS THE INTEGER REPRESENTATION
!                             USUALLY 32 BITS --> MAXINT=2147482648
!                             (ABOUT 68 YEARS).

!     METHOD.
!     -------

!       THE DATE AND TIME CAN BE SUPPLIED IN THE Y2K COMPLIANT
!       14 CHARACTER OR 12 CHARACTER FORMAT. IF THE 10 CHARACTER
!       FORMAT IS USED, AN ERROR MESSAGE WILL BE ISSUED.

!     EXTERNALS.
!     ----------

!     REFERENCES.
!     -----------


! ----------------------------------------------------------------------

      INTEGER :: MON(12)
      INTEGER :: IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC

      CHARACTER*(*) CDATE
      CHARACTER* 80 CLFMT
      INTEGER ISHIFT
      LOGICAL LLND

      DATA MON /31,28,31,30,31,30,31,31,30,31,30,31/

! ----------------------------------------------------------------------

!*    1.0 SPLIT DATE TIME GROUP INTO SECOND, MINUTE, HOUR, DAY, MONTH,
!*        YEAR DEPENDENT ON THE FORMAT OF CDATE.
!         -----------------------------------------------------------

      IL = LEN_TRIM(CDATE)
      IF (IL==10) THEN
        WRITE(CLFMT,'(A49, I2.2, A2)') &
     &     '(" INCDATE : NON-Y2K COMPLIANT  DATE IN USE ", A', IL, ' )'
        WRITE(6, CLFMT ) CDATE
        WRITE(*, CLFMT ) CDATE
        LLND = .FALSE.
        !CALL ABORT1
        STOP
      ELSEIF (IL==12) THEN
        LLND = .TRUE.
        READ(CDATE,'(i4, 4i2)')IYEAR, IMON, IDAY, IHOUR, IMIN
        ISEC = 0
      ELSEIF (IL==14) THEN
        LLND = .TRUE.
        READ(CDATE,'(i4, 5i2)')IYEAR, IMON, IDAY, IHOUR, IMIN, ISEC
      ELSEIF (IL==0) THEN
        WRITE(6, *) ' INCDATE:  FATAL@! DATE IS 0 CHARACTER LONG'
        WRITE(*, *) ' INCDATE:  FATAL@! DATE IS 0 CHARACTER LONG'
        !CALL ABORT1
        STOP
      ELSE
        WRITE(CLFMT,'(A59, I2.2, A2)') &
     &  '(" INCDATE:  FATAL@! DATE IS ",I2," CHARACTERS LONG!!= ",A', IL, ' )'
        WRITE(6, CLFMT) IL, CDATE
        WRITE(*, CLFMT) IL, CDATE
        STOP
      ENDIF

      IRET=0
      MON(2)=MFEB_LENGTH(IYEAR)
      IF(IMON<1 .OR. IMON>12) IRET=-1
      IF(IDAY<1 .OR. IDAY>MON(MIN(MAX(IMON,1),12))) IRET=-2
      IF(IHOUR<0 .OR. IHOUR>23) IRET=-3
      IF(IMIN<0 .OR. IMIN>59) IRET=-4
      IF(ISEC<0 .OR. ISEC>59) IRET=-5
      IF (IRET/=0) THEN
        WRITE(6, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(6, '("           INPUT DATE INCORRECT,")')
        WRITE(6, '("           IYEAR ", I4 )') IYEAR
        WRITE(6, '("           IMON ", I2 )') IMON
        WRITE(6, '("           IDAY ", I2 )') IDAY
        WRITE(6, '("           IHOUR ", I2 )') IHOUR
        WRITE(6, '("           IMIN ", I2 )') IMIN
        WRITE(6, '("           ISEC ", I2 )') ISEC
        WRITE(*, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(*, '("           INPUT DATE INCORRECT,")')
        WRITE(*, '("           IYEAR ", I4 )') IYEAR
        WRITE(*, '("           IMON ", I2 )') IMON
        WRITE(*, '("           IDAY ", I2 )') IDAY
        WRITE(*, '("           IHOUR ", I2 )') IHOUR
        WRITE(*, '("           IMIN ", I2 )') IMIN
        WRITE(*, '("           ISEC ", I2 )') ISEC
        !CALL ABORT1
        STOP
      ENDIF
      IF (IL==12 .AND. MOD(ISHIFT,60).NE.0) THEN
        IRET=-6
        WRITE(6, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(6, '(" WHEN HANDELING CHARACTER*12 DATES,")')
        WRITE(6, '(" ISHIFT MUST BE A MULTIPLE OF 60 SECONDS")')
        WRITE(6, '(" ISHIFT = ", I12 )') ISHIFT
        WRITE(*, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(*, '(" WHEN HANDELING CHARACTER*12 DATES,")')
        WRITE(*, '(" ISHIFT MUST BE A MULTIPLE OF 60 SECONDS")')
        WRITE(*, '(" ISHIFT = ", I12 )') ISHIFT
        !CALL ABORT1
        STOP
      ENDIF

! ----------------------------------------------------------------------

!*    2.0 INCREMENT OR DECREMENT THE DATE AND TIME BY ISHIFT SECONDS.
!         -----------------------------------------------------------


      ISEC=ISEC+MOD(ISHIFT,60)
      IF (ISEC.GE.60) THEN
        IMIN = IMIN+ISEC/60
        ISEC = ISEC-(ISEC/60)*60
      ELSE IF (ISEC.LT.0) THEN
        IMIN = IMIN +(ISEC-59)/60
        ISEC = ISEC-((ISEC-59)/60)*60
      END IF

      IMIN=IMIN+ISHIFT/60

!     2.1 POSITIVE SHIFT GREATER THAN 1 MINUTE.

      IF (IMIN.GE.60) THEN
         IHOUR = IHOUR + IMIN/60
         IMIN = IMIN - (IMIN/60)*60
         IF (IHOUR.GE.24) THEN
            IDAY = IDAY + IHOUR/24
            IHOUR = IHOUR - (IHOUR/24)*24
            DO WHILE (IDAY.GT.MON(IMON))
               IDAY=IDAY-MON(IMON)
               IMON=IMON+1
               IF(IMON.EQ.13) THEN
                  IMON = 1
                  IYEAR=IYEAR+1
                  MON(2)=MFEB_LENGTH(IYEAR)
               END IF
            ENDDO
         END IF
      ELSE IF (IMIN.LT.0) THEN

!     2.2 NEGATIVE SHIFT.

         IHOUR = IHOUR + (IMIN-59)/60
         IMIN = IMIN - ((IMIN-59)/60)*60
         IF (IHOUR.LT.0) THEN
            IDAY = IDAY + (IHOUR-23)/24
            IHOUR = IHOUR - ((IHOUR-23)/24)*24
            DO WHILE (IDAY.LT.1)
               IMON=IMON-1
               IF(IMON.EQ.0) THEN
                  IMON = 12
                  IYEAR=IYEAR-1
                  MON(2)=MFEB_LENGTH(IYEAR)
               END IF
               IDAY=IDAY+MON(IMON)
            ENDDO
         END IF
      END IF


! ----------------------------------------------------------------------

!*    3.0 CREATE OUTPUT STRING OF NEW DATE AND TIME USING THE INPUT
!*        FORMAT.
!         ---------------------------------------------------------

      IRET=0
      MON(2)=MFEB_LENGTH(IYEAR)
      IF(IMON<1 .OR. IMON>12) IRET=-1
      IF(IDAY<1 .OR. IDAY>MON(MIN(MAX(IMON,1),12))) IRET=-2
      IF(IHOUR<0 .OR. IHOUR>23) IRET=-3
      IF(IMIN<0 .OR. IMIN>59) IRET=-4
      IF(ISEC<0 .OR. ISEC>59) IRET=-5
      IF (IRET/=0) THEN
        WRITE(6, '(" INCDATE:  FATAL@! IRET= ", I3 )') IRET
        WRITE(6, '("           OUTPUT DATE INCORRECT,")')
        WRITE(6, '("           IYEAR ", I4 )') IYEAR
        WRITE(6, '("           IMON ", I2 )') IMON
        WRITE(6, '("           IDAY ", I2 )') IDAY
        WRITE(6, '("           IHOUR ", I2 )') IHOUR
        WRITE(6, '("           IMIN ", I2 )') IMIN
        WRITE(6, '("           ISEC ", I2 )') ISEC
        WRITE(6, '("           ISHIFT ", I12 )') ISHIFT
        !CALL ABORT1
        STOP
      ENDIF

      IF (LLND) THEN
        IF (IL==12) THEN
          WRITE(CDATE,'(I4.4, 4I2.2)')IYEAR, IMON, IDAY, IHOUR, IMIN
        ELSEIF (IL==14) THEN
          WRITE(CDATE,'(I4.4, 5I2.2)')IYEAR, IMON, IDAY, IHOUR,IMIN, ISEC
        ELSE
          WRITE(6,'(" THE LENGTH OF THE INPUT CHARACTER CHANGED @!")')
          !CALL ABORT1
          STOP
        ENDIF
      ELSE
        IYEAR = IYEAR - 1900
        IF (IYEAR >= 100) THEN
          WRITE(CLFMT,'(A39, I2.2, A2)') '(" NON-Y2K COMPLIANT  DATE IN USE ", A', IL, ' )'
          WRITE(6, CLFMT ) CDATE
          WRITE(6,'(" THIS CANNOT BE DONE ANY LONGER")')
          !CALL ABORT1
          STOP
        ENDIF
        WRITE(CDATE,'(5i2.2)')IYEAR, IMON, IDAY, IHOUR, IMIN
      ENDIF

      RETURN

      CONTAINS

      FUNCTION MFEB_LENGTH(IYEAR)
!     LENGTH OF FEBRUARY IN DAYS FOR YEAR IYEAR (YYYY)
      INTEGER :: MFEB_LENGTH
      INTEGER :: IYEAR
      IF(MOD(IYEAR,400).EQ.0) THEN
        MFEB_LENGTH=29
      ELSE IF(MOD(IYEAR,100).EQ.0) THEN
        MFEB_LENGTH=28
      ELSE IF(MOD(IYEAR,4).EQ.0) THEN
        MFEB_LENGTH=29
      ELSE
        MFEB_LENGTH=28
      ENDIF
      END FUNCTION MFEB_LENGTH

      END SUBROUTINE INCDATE

!!!!! FUNCTION julday !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Returns the Julian day number that begins at noon of the calendar date
! specified by month mm, day id, and year iyyy, all integers. Negative year
! signifies B.C.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION julday(iyyy,mm,id)

      IMPLICIT NONE

      ! Interface
      INTEGER iyyy, mm, id ! INTENT (IN)
      INTEGER julday       ! INTENT (OUT)

      ! Constant
      INTEGER IGREG
      PARAMETER ( IGREG=15+31*(10+12*1582) )

      ! Locals
      INTEGER ja, jm, jy

      ! Main
      jy = iyyy

      if (jy .EQ. 0) then
        print *, "JULDAY: There is no year zero"
        STOP
      endif

      if (jy .LT. 0) then
        jy = jy+1
      endif

      if (mm .GT. 2) then
        jm = mm+1
      else
        jy = jy-1
        jm = mm+13
      endif

      julday = 365*jy+int(0.25*jy+2000.0)+int(30.6001*jm)+id+1718995

      if (id+31*(mm+12*iyyy) .GE. IGREG) then
        ja = int(0.01*jy)
        julday = julday+2-ja+int(0.25*ja)
      endif

      END ! FUNCTION julday

!
! ---------------------------------------------------------------------------- !
