# RAD_AS_REF_4QDOAS

Python script for calculating averaged reference radiance from NADIR mode, band 2B for radiance as reference calculations in QDOAS.  Ideally use Python3 or higher.

Usage:

python3 rad_as_ref.py -i <inputfiles.nc> -o <outputfile.nc> or <outputfile.nc> --lat <start_lat> <end_lat> --lon <start_lon> <end_lon> (--apex or --text) (--irrad or --wlen <wavelength_start> <wavelength_end> --step <step>) --plot > output.txt

-i : Names of the input files. Supported input file format is NetCDF4 for GOME1. Multiple files are allowed, the script calculates the total average for all of them.

-o : Name of the output file where average is stored. Supported formats are APEX (extension is .nc) and ASCII (extension is .dat).

--lat : Latitudes within which the reference is calculated.

--lon : Longitudes within which reference is calculated.

--apex or --text : formats of the output files.

--irrad : Irradiance spectral channel of the first relevant input file chosen as the interpolation grid. Not compatible with following two options.

--wlen : Wavelength range of the interpolation grid. Not compatible with --irrad.

--step : Interpolation step in nanometers for the interpolation grid. Necessary if --wlen defined and not compatible with --irrad.

--plot : Plot the reference spectra.

output.txt: name of the file for verbose printing of the workflow sent to the standard output.
