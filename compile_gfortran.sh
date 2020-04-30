#!/bin/bash

FLAGS="-O2 -mavx2 -m64 -march=native"
# FLAGS="$FLAGS -pg"

if [ ! -d build ]; then
    mkdir build
fi

cd build/
pwd
rm -fv *.o *.mod

export GRIB_API_INCLUDE="-I/usr/local/include"
export GRIB_API_LIB="-lgrib_api_f90 -lgrib_api"

export NETCDF_INCLUDE="$(nf-config --fflags)"
export NETCDF_LIB="$(nf-config --flibs) -lnetcdff"

rm -fv *.mod
rm -fv *.o

gfortran $FLAGS $NETCDF_INCLUDE $NETCDF_LIB -lnetcdff -c ../src/netcdf_metno_spec.f90
gfortran $FLAGS -c ../src/debug.f90
gfortran $FLAGS -c ../src/sphere_vec.f90
gfortran $FLAGS -fopenmp ../src/debug.f90 ../src/gribspec2nc.f90 ../src/sphere.f ../src/sphere_vec.f90 ../src/netcdf_metno_spec.f90 $NETCDF_INCLUDE $NETCDF_LIB $GRIB_API_INCLUDE $GRIB_API_LIB -o ../gribspec2nc

cd ..
