#!/bira-iasb/softs/18/py36/bin/python3.6

import os
import sys
import math
import datetime as dt
import time
import numpy as np
import numpy.ma as ma
import pandas as pa
import string
import argparse as ap
#from netCDF4 import Dataset
import netCDF4
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import ascii
from io import StringIO
from GOME_RadAsRef import mainGOME
from TROPOMI_RadAsRef import mainTROPOMI

parser = ap.ArgumentParser(description='Script that calculates a daily mean of radiances within given coordinates for the RadAsRef option in QDOAS',
                           epilog='''Thank you for using it! 
                                     Dalibor Hršak''',
                           usage='%(prog)s -i <inputfiles.nc> -o <outputfile.nc> or <outputfile.dat> --lat <start_lat> <end_lat> --lon <start_lon> <end_lon> --irrad or --wlen <wlen_start> <wlen_end> <step> (--apex or --text) --plot',
                           fromfile_prefix_chars='@')

parser.add_argument('-i', dest='inputfiles', nargs='+', metavar='INPUT_FILE',
                    default=[],
                    help='''Specify the names of the input NetCDF4 files''')

parser.add_argument('-o', dest='outputfiles', nargs=1, metavar='OUTPUT_FILE',
                    default=[],
                    help='''Specify the name of the output file.
                            For ASCII format extension is .dat, for NetCDF4 it is .nc.''')

parser.add_argument('--gome', action='store_true',
                    help='''Apply the program on a set of GOME data.''')

parser.add_argument('--tropo', action='store_true',
                    help='''Apply the program on a set of TROPOMI data.''')

parser.add_argument('--lat',nargs=2,default=[],dest='latitude',type=float,
                    metavar=('LAT_START','LAT_END'),
                    help='''Give latitude range from -90 to 90.''')

parser.add_argument('--lon',nargs=2,default=[],dest='longitude',type=float,
                    metavar=('LON_START','LON_END'),
                    help='''Give longitude range from 0 to 360.''')

parser.add_argument('--wlen',nargs=3,default=[],dest='wavelength',type=float,
                    metavar=('WLEN_START','WLEN_END','STEP'),
                    help='''Give wavelength range of interest in nm for the custom interpolation grid and the step.''')

parser.add_argument('--irrad', action='store_true',
                    help='''Choose the irradiance wavelengths as the interpolation grid.''')

parser.add_argument('--apex', action='store_true',
                    help='''Write the reference spectra to a netCDF4 file in the APEX format.''')

parser.add_argument('--text', action='store_true',
                    help='''Write the reference spectra to a TEXT file in the ASCII format.''')

parser.add_argument('--plot', action='store_true',
                    help='''Plot the reference spectra in a PDF format.''')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit()
args = parser.parse_args()
do_plot=args.plot
do_irrad=args.irrad
do_ascii=args.text
do_apex=args.apex
gome=args.gome
tropomi=args.tropo

if not (args.inputfiles and args.outputfiles and (args.gome or args.tropo) and args.latitude and args.longitude):
    parser.print_help()
    sys.exit()

if (len(args.latitude) != 2 or len(args.longitude) != 2):
    print('ERROR: latitude and longitude need exactly 2 arguments!')
    sys.exit()

if (len(args.wavelength) != 3 and not do_irrad and not tropomi):
    print('ERROR: keyword wavelength needs exactly 2 arguments!')
    sys.exit()

if (len(args.latitude) < 2 or len(args.longitude) < 2):
    print('ERROR: too few arguments for latitude or longitude (min. 2)')
    sys.exit()

if (len(args.latitude) > 2 or len(args.longitude) > 2):
    print('ERROR: too many arguments for latitude or longitude (max. 2)')
    sys.exit()

if (gome and tropomi):
    print('ERROR: either choose --gome or --tropo!')
    sys.exit()

if (gome and not(do_apex or do_ascii)):
    print('ERROR: For GOME choose APEX or ASCII format!')
    sys.exit()

if (gome and not(do_irrad or args.wavelength)):
    print('ERROR: For GOME choose irradiance wavelength interpolation grid or define your own!')
    sys.exit()

if (args.wavelength and do_irrad):
    print('ERROR: either define the interpolation grid or choose the irradiance grid!')
    sys.exit()

if (do_ascii and do_apex):
    print('ERROR: either choose --apex or --text!')
    sys.exit()

if (do_ascii and tropomi):
    print('''WARNING: Storing TROPOMI reference in ASCII format not available!
             Proceeding, reference will be stored in NetCDF4 format.''')

if (tropomi and do_irrad):
    print('''WARNING: For TROPOMI do irradiance not needed!
             Proceeding without interpolation.''')

if (do_plot and tropomi):
    print('''WARNING: plotting not available for TROPOMI!
             Proceeding without plotting.''')
    do_plot = False

longitude = args.longitude
latitude = args.latitude
if (args.longitude[0] < args.longitude[1]):
    longs = 1
elif (args.longitude[0] > args.longitude[1]):
    longs = 2
for i in range(len(args.longitude)):
    if (args.longitude[i] > 180.0 and args.longitude[i] < 360.0):
        longitude[i] = args.longitude[i] - 360.0
    elif (args.longitude[i] >= 0.0 and args.longitude[i] <= 180.0):
        longitude[i] = args.longitude[i]
    else:
        print('ERROR: longitude outside of range (0 to 360 degrees)!')
        sys.exit()

start_total = time.time()
if (gome):
    if (do_irrad):
        mainGOME(args.inputfiles,args.outputfiles[0],latitude,longitude,longs,do_irrad,do_apex,do_plot)
    else:
        mainGOME(args.inputfiles,args.outputfiles[0],latitude,longitude,longs,do_irrad,do_apex,do_plot,args.wavelength)
elif (tropomi):
    if (args.wavelength):
        mainTROPOMI(args.inputfiles,args.outputfiles[0],latitude,longitude,longs,do_apex,args.wavelength)
    else:
        mainTROPOMI(args.inputfiles,args.outputfiles[0],latitude,longitude,longs,do_apex)
end_total = time.time()
print('>---Total elapsed time for today = {0:7.5f} seconds.---<\n'.format(end_total-start_total))
print('%'*(36)+'\n'+'%% Closing the program for today. %%\n'.format(args.outputfiles[0])+'%'*(36))
