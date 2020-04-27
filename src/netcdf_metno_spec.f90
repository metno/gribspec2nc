!########################################################
!# Module for creating NetCDF file with 2D wave spectra
!# and integrated parameters from selected locations
!#
!# Oyvind Breivik
!# oyvind.breivik@met.no
!# 2016-01-14
!########################################################
MODULE netcdf_metno_spec

   USE netcdf

   IMPLICIT NONE
   PRIVATE

   INTEGER, PARAMETER :: NFS=6             !! maximum number of fields, added 3 scalar fields
   INTEGER, PARAMETER :: IU06=6

   INTEGER, SAVE :: ncidspc(0:NFS+7) = -1 ! had to change from 3 to 7.

   PUBLIC nc_open_specfile, nc_write_spec, nc_close_specfile, nc_write_intpar

   CONTAINS 


   !########################################################
   !#      write header of a spectral NetCDF file    
   !########################################################
   SUBROUTINE nc_open_specfile(fname, sd, cflag_s, speclong, speclat, fr, th, ideldo)
    
      IMPLICIT NONE

      !CHARACTER (LEN=16), INTENT(IN) :: fname
      CHARACTER (LEN=128), INTENT(IN) :: fname
      CHARACTER (LEN=14), INTENT(IN) :: sd
      LOGICAL,  INTENT(IN) :: cflag_s(:)
      REAL,     INTENT(IN) :: speclong(:)  
      REAL,     INTENT(IN) :: speclat(:)
      REAL,     INTENT(IN) :: fr(:)
      REAL,     INTENT(IN) :: th(:)
      INTEGER,  INTENT(IN) :: ideldo
    
      CHARACTER (LEN=33) :: tua
      CHARACTER (LEN=21) :: tda
      CHARACTER (LEN=42) :: metno
      CHARACTER (LEN=100), DIMENSION (4, nfs) :: vl  
      INTEGER, DIMENSION (5)   :: disd
      INTEGER, DIMENSION (0:4) :: dsd
      INTEGER, DIMENSION (3)   :: diid
      INTEGER, DIMENSION (0:2) :: did
      INTEGER, DIMENSION (nfs)  :: flg, ty
      INTEGER :: n, m, l, i, loop, ii
      REAL*8 :: dtor, dt
      REAL*4, DIMENSION (nfs) :: fw, vmin, vmax
      INTEGER :: yy
      INTEGER, ALLOCATABLE, DIMENSION (:) :: xx
      
      metno = 'Norwegian Meteorological Institute, met.no'
      ALLOCATE (xx(size(speclong)))
      xx = (/ (ii, ii=1,size(speclong)) /)
      yy = 1
        
      dt = real(ideldo)
      IF (ncidspc(0)<0) THEN
         vl  = ''                        
         vl(1, 1) = 'SPEC'
         vl(2, 1) = '2D_SPECTRUM'
         vl(3, 1) = '2-D spectrum of total sea'
         vl(4, 1) = 'm**2 s  '
         vmin(1)  =  0.00
         vmax(1)  =  100.0
        
         vl(1, 2) = 'SPEC_SEA'
         vl(2, 2) = '2D_SPECTRUM_OF_WIND_SEA'
         vl(3, 2) = 'Wind Sea 2-D spectrum'
         vl(4, 2) = 'm**2 s'
         vmin(2)  =  0.0
         vmax(2)  =  100.0
    
         vl(1, 3) = 'SPEC_SWELL'
         vl(2, 3) = '2D_SPECTRUM_OF_WIND_SWELL'
         vl(3, 3) = 'Swell 2-D spectrum'
         vl(4, 3) = 'm**2 s'
         vmin(3)  =  0.0
         vmax(3)  =  100.0

         vl(1, 4) = 'hs'
         vl(2, 4) = 'sea_surface_wave_significant_height'
         vl(3, 4) = 'Total significant wave height'
         vl(4, 4) = 'm'
         vmin(4)  = 0.
         vmax(4)  = 20.

         vl(1, 5) = 'tp'
         vl(2, 5) = 'sea_surface_wave_peak_period_from_variance_spectral_density'
         vl(3, 5) = 'Total peak period'
         vl(4, 5) = 's'
         vmin(5) = 1.
         vmax(5) = 30.
      
         vl(1, 6) = 'Pdir'
         vl(2, 6) = 'peak_direction'
         vl(3, 6) = 'peak direction'
         vl(4, 6) = 'degree'
         vmin(6) = 0.
         vmax(6) = 360.

         flg = 1                   !! prepare table for required parameters only
         DO loop = 1,nfs
            IF (.NOT. cflag_s(loop)) THEN
               flg([loop]) = 0
            ENDIF
         ENDDO
      
         ty = NF90_FLOAT           !! type of parameter
         fw = -99999.0             !! dummy (zmiss)
         fw(4:6) = -999.0          !! dummy (zmiss)
         dsd = [5,4,3,2,1]
         did = [3,2,1]
         i = dt
         WRITE (tda,'("0000-00-00 (",2(i2.2,":"),i2.2,")")') i/3600,MOD(i/60,60),MOD(i,60)
         tua="seconds since "//sd(1:4)//"-"//sd(5:6)//"-"//sd(7:8)//" "//sd(9:10)//":"//sd(11:12)//":"//sd(13:14)
         CALL pf(nf90_create(fname, ior(NF90_CLOBBER,NF90_SHARE), ncidspc(0)))
         CALL pf(nf90_def_dim(ncidspc(0), 'time', NF90_UNLIMITED, disd(1))) 
         CALL pf(nf90_def_dim(ncidspc(0), 'x', size(xx) ,disd(3)))
         CALL pf(nf90_def_dim(ncidspc(0), 'y', 1, disd(2)))
         CALL pf(nf90_def_dim(ncidspc(0), 'freq', size(fr), disd(4)))
         CALL pf(nf90_def_dim(ncidspc(0), 'direction', size(th), disd(5)))
         CALL pf(nf90_def_var(ncidspc(0), 'direction', NF90_FLOAT, [disd(5)],ncidspc(nfs+1)))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+1), 'units', 'degree'))
         CALL pf(nf90_def_var(ncidspc(0), 'freq', NF90_FLOAT, [disd(4)], ncidspc(nfs+2)))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+2), 'units', '1/s'))
         CALL pf(nf90_def_var(ncidspc(0), 'x', NF90_INT, [disd(3)], ncidspc(nfs+3)))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+3), 'axis', 'x'))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+3), 'long_name', 'x-coordinate in Cartesian system'))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+3), 'standard_name', 'projection_x_coordinate')) 
         CALL pf(nf90_put_att(ncidspc(0),ncidspc(nfs+3),'units','1'))
         CALL pf(nf90_def_var(ncidspc(0),'y',nf90_int,[disd(2)], ncidspc(nfs+4)))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+4), 'axis', 'y'))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+4), 'long_name', 'y-coordinate in Cartesian system'))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+4), 'standard_name','projection_y_coordinate')) 
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+4), 'units','1'))
         CALL pf(nf90_def_var(ncidspc(0), 'time',NF90_INT, [disd(1)],ncidspc(nfs+5)))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+5), 'delta_t', tda))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+5), 'units',tua))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+5), 'dt', i))
         CALL pf(nf90_def_var(ncidspc(0), 'latitude', nf90_float, [disd(3), disd(2)],ncidspc(nfs+6)))  
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+6), 'long_name', 'latitude'))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+6), 'standard_name', 'latitude'))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+6), 'units', 'degree_north'))
         CALL pf(nf90_def_var(ncidspc(0), 'longitude', NF90_FLOAT, [disd(3), disd(2)],ncidspc(nfs+7)))  
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+7), 'long_name', 'longitude'))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+7), 'standard_name','longitude'))
         CALL pf(nf90_put_att(ncidspc(0), ncidspc(nfs+7), 'units', 'degree_east'))
         CALL pf(nf90_put_att(ncidspc(0), nf90_global, 'Conventions', 'CF-1.0'))
         CALL pf(nf90_put_att(ncidspc(0), nf90_global, 'institution', metno))
         
         ! Spectra
         DO i = 1,3
            IF (flg(i)>0) THEN
               CALL pf(nf90_def_var(ncidspc(0), vl(1,i), ty(i), disd(dsd), ncidspc(i)))
               CALL pf(nf90_put_att(ncidspc(0), ncidspc(i), '_FillValue', fw(i)))
               CALL pf(nf90_put_att(ncidspc(0), ncidspc(i), 'long_name', vl(3,i)))
               CALL pf(nf90_put_att(ncidspc(0), ncidspc(i), 'standard_name', vl(2,i)))
               CALL pf(nf90_put_att(ncidspc(0), ncidspc(i), 'units', vl(4,i)))
            ENDIF
         ENDDO

         ! Integrated parameters
         DO i = 4, 6
            IF (flg(i)>0) THEN
               CALL pf(nf90_def_var(ncidspc(0), vl(1,i), ty(i), disd(did), ncidspc(i)))
               CALL pf(nf90_put_att(ncidspc(0), ncidspc(i), '_FillValue', fw(i)))
               CALL pf(nf90_put_att(ncidspc(0), ncidspc(i), 'long_name', vl(3,i)))
               CALL pf(nf90_put_att(ncidspc(0), ncidspc(i), 'standard_name', vl(2,i)))
               CALL pf(nf90_put_att(ncidspc(0), ncidspc(i), 'units', vl(4,i)))
               CALL pf(nf90_put_att(ncidspc(0), ncidspc(i), 'valid_min', vmin(i)))
               CALL pf(nf90_put_att(ncidspc(0), ncidspc(i), 'valid_max', vmax(i)))
            ENDIF
         ENDDO
        
         CALL pf(nf90_enddef(ncidspc(0)))
         CALL pf(nf90_put_var(ncidspc(0),ncidspc(nfs+1),th))
         CALL pf(nf90_put_var(ncidspc(0),ncidspc(nfs+2),fr))
         CALL pf(nf90_put_var(ncidspc(0),ncidspc(nfs+3),xx))
         CALL pf(nf90_put_var(ncidspc(0),ncidspc(nfs+4),yy))
         CALL pf(nf90_put_var(ncidspc(0),ncidspc(nfs+6),speclat))
         CALL pf(nf90_put_var(ncidspc(0),ncidspc(nfs+7),speclong))
      ELSE
         write (IU06,*) ' +++ NetCDF-file is open already'
         STOP  ' +++ Problem occurs during NetCDF-Output handling '
      ENDIF
       
   END SUBROUTINE nc_open_specfile 
   ! ---------------------------------------------------------------------------- !
      
   !########################################################
   !# 	Write spectra to NetCDF-file
   !########################################################
   SUBROUTINE nc_write_spec(ip, spec, t, ipos, tda)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ip          !! parameter number
      REAL,    INTENT(IN) :: spec(:,:)   !! spectrum
      INTEGER, INTENT(IN) :: t
      INTEGER, INTENT(IN) :: ipos    
      INTEGER, INTENT(IN) :: tda

      INTEGER, DIMENSION (5) :: start, icnt
      INTEGER :: vid, j !, tda

      INTEGER, SAVE :: no = -1                !! offset
      LOGICAL, SAVE :: fc = .TRUE.            !! first call?
      
      IF (fc) THEN
         no = 1-t                           !! t+no=1 : first output
         fc = .FALSE.
      ENDIF
      vid = ip
      
      IF (ncidspc(vid)>0) THEN
         start = [1,   1, ipos , 1, t+no]
         icnt   = [size(spec,1), size(spec,2), 1 , 1, 1]
         CALL pf(NF90_PUT_VAR(ncidspc(0),ncidspc(vid),spec,start,icnt))
         CALL pf(NF90_PUT_VAR(ncidspc(0),ncidspc(nfs+5),t*tda,[t+no]))
      ELSE
         WRITE (IU06,*) ' +++ parameter ',vid, 'is not available in NetCDF-file'
      ENDIF
      
   END SUBROUTINE nc_write_spec


   !########################################################
   !# 	Write additional integrated parameters to NetCDF-file
   !########################################################
   SUBROUTINE nc_write_intpar(vid, grid, t, tda)
    
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: vid         !! parameter number
      REAL,    INTENT(IN) :: grid(:,:)   !! grid of parameter
      INTEGER, INTENT(IN) :: t
      INTEGER, INTENT(IN) :: tda

      INTEGER, DIMENSION (3) :: sta
      INTEGER :: j

      INTEGER, SAVE :: no = -1                !! offset
      LOGICAL, SAVE :: fc = .TRUE.            !! first call?
      
      IF (fc) THEN
         no = 1-t                           !! t+no=1 : first output
         fc = .FALSE.
      ENDIF
      
      IF (ncidspc(vid)>0) THEN
         sta = [1,1,t+no]
         CALL pf(NF90_PUT_VAR(ncidspc(0),ncidspc(vid),grid,sta))
         CALL pf(NF90_PUT_VAR(ncidspc(0),ncidspc(nfs+5),t*tda,[t+no]))
      ELSE
         WRITE (IU06,*) ' +++ parameter ',vid, 'is not available in NetCDF-file'
      ENDIF
      
   END SUBROUTINE nc_write_intpar


   !##
   !##############################
   !# 	close NetCDF-file
   !##############################
   SUBROUTINE nc_close_specfile
      IF (ncidspc(0)>=0) CALL pf(NF90_CLOSE(ncidspc(0)))
   END SUBROUTINE nc_close_specfile

      
   !################################################
   !# pf write a NetCDF-error message 
   !#	en	NetCDF error number         
   !################################################
   SUBROUTINE pf(en)
      INTEGER :: en
      IF (en/=0) WRITE (IU06,*) ' +++ Error : ', NF90_STRERROR(en)
   END SUBROUTINE pf

END MODULE netcdf_metno_spec
