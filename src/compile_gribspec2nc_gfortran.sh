#!/bin/bash

export GRIB_API_INCLUDE="-I/usr/local/include"
export GRIB_API_LIB="-L/usr/local/lib -lgrib_api_f90 -lgrib_api"

export NETCDF_INCLUDE="$(nf-config --fflags)"
export NETCDF_LIB="$(nf-config --flibs) -lnetcdff"

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
echo $LD_LIBRARY_PATH

rm -fv *.mod
rm -fv *.o

gfortran -c debug.f90
gfortran $NETCDF_INCLUDE $NETCDF_LIB -lnetcdff -c netcdf_metno_spec.f90
gfortran sphere.f debug.f90 gribspec2nc.f90 netcdf_metno_spec.f90 $NETCDF_INCLUDE $NETCDF_LIB $GRIB_API_INCLUDE $GRIB_API_LIB -o gribspec2nc
