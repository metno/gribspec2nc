# GribSpec2NC

Oyvind Breivik, Norwegian Meteorological Institute

## PURPOSE

Decodes GRIB spectra stored as GRIB parameter 251.
The results will be saved as a NetCDF file containing
spectra in selected locations provided as input.
The spectra are interpolated to these locations if possible.

## USAGE

```
    gribspec2nc [-x] [-d timestep_hh][-i infile] [-l speclist.inp] [-o outfile] [-t itest] [-w weightfile]

    -x: compute integrated parameters
    -d timestep_hh: specify temporal resolution, default is 1 h
    -i infile: GRIB input, default is input_spectra
    -l speclist.inp: list of spectral locations to interpolate to,
                default is speclist.inp
    -o outfile: NetCDF output, default is output_spectra.nc
    -w weightfile: Use precomputed interpolation weights, default is no
    -s: issue a warning but generate NetCDF file even if there are missing GRIB fields (default is to issue error and abort)
    -t: debug reporting:
        itest: 0 (no diagnostics)
        itest: 1 (some diagnostics)
        itest: 2 (some more diagnostics)
        itest: 3 (all diagnostics)
```

## INPUT FILE REQUIREMENT

The input file can only contain one parameter (251) obtained
for all directions and frequencies in the default order. Namely, the
mars request should be done with DIRECTION=1/TO/nang, and
FREQUENCY=1/TO/nfre. At this time nfre=30 and nang=24 for global model
and for mediterranean data.
