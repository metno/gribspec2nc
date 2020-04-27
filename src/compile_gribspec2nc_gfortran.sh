#!/bin/bash

# FLAGS="-pg"
FLAGS="-O2"

export GRIB_API_INCLUDE="-I/usr/local/include"
export GRIB_API_LIB="-L/usr/local/lib -lgrib_api_f90 -lgrib_api"

export NETCDF_INCLUDE="$(nf-config --fflags)"
export NETCDF_LIB="$(nf-config --flibs) -lnetcdff"

rm -fv *.mod
rm -fv *.o

gfortran $FLAGS $NETCDF_INCLUDE $NETCDF_LIB -lnetcdff -c netcdf_metno_spec.f90
gfortran $FLAGS -c sphere_vec.f90
gfortran $FLAGS gribspec2nc.f90 sphere.f sphere_vec.f90 netcdf_metno_spec.f90 $NETCDF_INCLUDE $NETCDF_LIB $GRIB_API_INCLUDE $GRIB_API_LIB -o gribspec2nc
# gfortran $FLAGS gribspec2nc.f90 netcdf_metno_spec.f90 $NETCDF_INCLUDE $NETCDF_LIB $GRIB_API_INCLUDE $GRIB_API_LIB -o gribspec2nc
