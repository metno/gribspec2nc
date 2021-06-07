#!/bin/bash

FLAGS="-O2 -mavx2 -m64 -march=native"
# FLAGS="$FLAGS -pg"

if [ ! -d build ]; then
    mkdir build
fi

cd build/
pwd
rm -fv *.o *.mod

# Load Modules
module purge
module load grib_api/1.28.0
module load hdf5/1.10.5-gcc
module load netcdf/4.7.0-gcc

#export GRIB_API_INCLUDE="-I/usr/local/include"
#export GRIB_API_LIB="-lgrib_api_f90 -lgrib_api"

export NETCDF_INCLUDE="$(nf-config --fflags)"
export NETCDF_LIB="$(nf-config --flibs) -lnetcdff"

gfortran $FLAGS $CPPFLAGS -c ../src/netcdf_metno_spec.f90
gfortran $FLAGS -c ../src/debug.f90 ../src/sphere_vec.f90
gfortran $FLAGS $CPPFLAGS $LDFLAGS ../src/debug.f90 ../src/gribspec2nc.f90 ../src/sphere.f ../src/sphere_vec.f90 ../src/netcdf_metno_spec.f90 -fopenmp -lnetcdff -lgrib_api_f90 -lgrib_api -o ../gribspec2nc

cd ..
