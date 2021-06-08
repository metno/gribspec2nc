#!/bin/bash

FLAGS="-O2 -mavx2 -m64 -march=native -fopenmp" # -std=f95"
# FLAGS="$FLAGS -pg" # Debug

if [ ! -d build ]; then
    mkdir build
fi

cd build/
pwd
rm -fv *.o *.mod

export FFLAGS="$(nf-config --fflags) -I/usr/lib/x86_64-linux-gnu/fortran/gfortran-mod-15"
export LDFLAGS="$(nf-config --flibs)"

gfortran $FLAGS $FFLAGS -c ../src/netcdf_metno_spec.f90
gfortran $FLAGS -c ../src/debug.f90 ../src/sphere_vec.f90 ../src/sphere.f
gfortran $FLAGS $FFLAGS $LDFLAGS ../src/gribspec2nc.f90 \
    debug.o sphere.o sphere_vec.o netcdf_metno_spec.o \
    -lnetcdff -leccodes_f90 -o ../gribspec2nc

cd ..
