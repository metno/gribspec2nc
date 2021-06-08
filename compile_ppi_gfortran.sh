#!/bin/bash

FLAGS="-O2 -mavx2 -m64 -march=native -fopenmp"
# FLAGS="$FLAGS -pg" # Debug

if [ ! -d build ]; then
    mkdir build
fi

cd build/
pwd
rm -fv *.o *.mod

# Load Modules
module purge
module load hdf5/1.10.5-gcc
module load netcdf/4.7.0-gcc
module load eccodes/2.20.0

gfortran $FLAGS $CPPFLAGS -c ../src/netcdf_metno_spec.f90
gfortran $FLAGS -c ../src/debug.f90 ../src/sphere_vec.f90
gfortran $FLAGS $CPPFLAGS $CFLAGS $LDFLAGS ../src/gribspec2nc.f90 \
    debug.o sphere_vec.o netcdf_metno_spec.o \
    -lnetcdff -leccodes_f90 -o ../gribspec2nc

cd ..
