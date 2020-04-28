#!/bin/bash

# purge modules to avoid conflicts
module purge

# load modules
module load buildenv-intel/2018a-eb
module load netCDF/4.4.1.1-HDF5-1.8.19-nsc1-intel-2018a-eb
module load grib_api/1.24.0-nsc1-intel-2018a-eb

# setup for grib
export GRIB_API=${GRIB_API_DIR}
export GRIB_API_INCLUDE="-I${GRIB_API}/include"
export GRIB_API_LIB="-L${GRIB_API}/lib -lgrib_api_f90 -lgrib_api"
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${GRIB_API}/lib
echo $LD_LIBRARY_PATH

# setup for netcdf
export NETCDF=${NETCDF_DIR}
export NETCDF_INCLUDE="-I${NETCDF}/include"
export NETCDF_LIB="-L${NETCDF}/lib -lnetcdff"

# compile module
ifort $NETCDF_INCLUDE $NETCDF_LIB -c netcdf_metno_spec.f90
# compile and generate executables
ifort gribspec2nc.f90 sphere.f netcdf_metno_spec.f90 $NETCDF_INCLUDE $NETCDF_LIB $GRIB_API_INCLUDE $GRIB_API_LIB -o gribspec2nc
