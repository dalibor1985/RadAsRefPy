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
from netCDF4 import Dataset
#import netCDF4
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import ascii
from io import StringIO
#from tabulate import tabulate
from ReadNC import GOMEData
from DataAnalysis import AnalyseGOME

def write_text(wavelengths,results,outfile):
    if (not outfile.endswith('.dat')):
        outfile = os.path.splitext(outfile)[0] + ".dat"
    output=''
    for i in range(len(wavelengths)):
        output += '{0:10.5f}{1:20.3f}{2:20.3f}{3:20.3f}{4:20.3f}\n'.format(wavelengths[i],results[0][i],results[1][i],results[2][i],results[3][i])
    wfile = open(outfile,'w')
    wfile.write(output)
    wfile.close()

def write_apex(wavelengths,results,outfile,numSpectra,lat,lon,inputfiles,fillValue):
    if (not outfile.endswith('.nc')):
        outfile = os.path.splitext(outfile)[0] + ".nc"
    longitude=np.empty(2)
    for i in range(len(lon)):
        if (lon[i] > -180.0 and lon[i] < 0.0):
            longitude[i] = lon[i] + 360.0
        elif (lon[i] >= 0.0 and lon[i] <= 180.0):
            longitude[i] = lon[i]
    year=inputfiles[0][20:24]
    month=inputfiles[0][24:26]
    day=inputfiles[0][26:28]
    dataset = Dataset(outfile, 'w',  format='NETCDF4_CLASSIC')
    dataset.title = 'Reference file created by RAD_AS_REF for QDOAS Python script.'
    dataset.description = 'Radiance as reference file in APEX format for QDOAS based on daily averaged radiances.'
    dataset.created = time.ctime(time.time())
    dataset.num_spectra = "{0} East, {1} Center, {2} West, {3} Backscan spectra.".format(numSpectra[0],numSpectra[1],numSpectra[2],numSpectra[3])
    dataset.lat_bound = lat[0],lat[1]
    dataset.lon_bound = longitude[0],longitude[1]
    dataset.inputfiles = inputfiles
    dataset.measurement_date = "{0}/{1}/{2}".format(day,month,year)
    dataset.fillValue = fillValue
    col_dim = dataset.createDimension('col_dim', np.shape(results)[0])
    spectral_dim = dataset.createDimension('spectral_dim', np.shape(wavelengths)[0])
    reference_radiance = dataset.createVariable('reference_radiance', np.float64, ('col_dim','spectral_dim'))
    reference_wavelength = dataset.createVariable('reference_wavelength', np.float64, ('spectral_dim',))
    use_row = dataset.createVariable('use_row', np.int8, ('col_dim',))
    reference_radiance[:] = results
    reference_wavelength[:] = wavelengths
    use_row[:] = np.where(numSpectra>0,1,0)
    dataset.close()

def plot_result(wavelengths,dataEast,dataCenter,dataWest,dataBack,outfile):
    plt.plot(wavelengths,dataEast,'-', linewidth=1, label='daily mean East') # plot the data points
    plt.plot(wavelengths,dataCenter,'-', linewidth=1, label='daily mean Center') # plot the data points
    plt.plot(wavelengths,dataWest,'-', linewidth=1, label='daily mean West') # plot the data points
    plt.plot(wavelengths,dataBack,'-', linewidth=1, label='daily mean backscan') # plot the data points
    plt.xlabel('wavelength (nm)')
    plt.ylabel('photons / (nm s cm^2)')
    plt.grid()
    plt.legend(loc='best')
    plt.title('Averaged spectrum')
    plotfile = os.path.splitext(outfile)[0] + '.pdf'
    plt.savefig(plotfile)

def transformLongitudes(inputLongitudes,inputMask):
    assert(np.shape(inputLongitudes) == np.shape(inputMask))
    longitudes = ma.array(inputLongitudes,mask=inputMask)
    for i in range(len(longitudes)):
        for j in range(len(longitudes[i])):
            if (longitudes[i][j] < 0.0 and longitudes[i][j] > -180.):
                longitudes[i][j] = longitudes[i][j] + 360.0
    return longitudes

def mainGOME(inputfiles,outputfile,latitude,longitude,longs,do_irrad,do_apex,do_plot,*wavelength):
#    assert(((do_irrad==True) and (len(wavelength[0]) == 0)) or ((do_irrad==False) and ((len(wavelength[1]) == 3))))
    # start the main loop over NetCDF4 files
    masterEast=[]
    masterCenter=[]
    masterWest=[]
    masterBack=[]

    print('%'*(44))
    print('%% RadAsRef.py script GOME mode for QDOAS. %%')
    print('%'*(44)+'\n')

    idx=0
    numSpectra=np.full(4,0)
    usedfiles = []
    for File in range(len(inputfiles)):
        start = time.time()
        filename=inputfiles[File]
        if (not filename.endswith('.nc')):
            filename += '.nc'
        if (not os.path.isfile(filename)):
            print('+'+'='*(30 + len(filename))+'+')
            print('| ERROR: file {0} does not exist! |'.format(filename))
            print('+'+'='*(30 + len(filename))+'+\n')
            sys.exit()
        else:
            print('+'+'='*(23 + len(filename))+'+')
            print('| Opening NetCDF4 file {0} |'.format(filename))
            print('+'+'='*(23 + len(filename))+'+\n')
        inputData = GOMEData()
        inputData.read(filename)
        end = time.time()
        # If MODE_NADIR is empty, stop and jump to the next input file
        if (inputData.empty_nadir):
            print('>---NADIR group in NetCDF4 file {0} is empty.---<'.format(filename))
            continue
        print('>---Time spent to read the NetCDF4 file = {0:7.5f} seconds.---<'.format(end-start))
        
        # Call the function that splits the data into East/Center/West/Backscan sets
        startSplit = time.time()
        splitData=AnalyseGOME()
        
        maskForward = splitData.maskGeodata(inputData.latitudes[0],inputData.longitudes[0],latitude,longitude,longs)
        maskBackward = splitData.maskGeodata(inputData.back_latitudes[0],inputData.back_longitudes[0],latitude,longitude,longs)
        latitudes = ma.array(inputData.latitudes[0],mask=maskForward)
        longitudes = transformLongitudes(inputData.longitudes[0],maskForward)
        back_latitudes = ma.array(inputData.back_latitudes[0],mask=maskBackward)
        back_longitudes = transformLongitudes(inputData.back_longitudes[0],maskBackward)
        #    np.set_printoptions(threshold=np.inf,precision=5)
        
        # if all data are masked, i.e. all coordinates are outside of desired range, stop and jump to the next input file
        if (np.array(maskForward).all() and np.array(maskBackward).all()):
            print('>---All masks are true, i.e. no points within given coordinates.---<')
            print('>---Continuing to the next file or exiting the loop.---<\n')
            continue
        print('>---Initially we have {0} spectra in forward scan.---<'.format(np.shape(inputData.radiance)[1]))
        print('>---Initially we have {0} spectra in back scan.---<'.format(np.shape(inputData.back_radiance)[1]))
        print('>---Now we split spectra into East, Center, West and back scan.---<')
        splitData.splitSpectra(inputData.radiance[0],maskForward)
        splitData.splitSpectra(inputData.back_radiance[0],maskBackward)
        fillValue=inputData.FillValue
        numSpectra[0]+=np.shape(splitData.dataEast)[0]
        numSpectra[1]+=np.shape(splitData.dataCenter)[0]
        numSpectra[2]+=np.shape(splitData.dataWest)[0]
        numSpectra[3]+=np.shape(splitData.dataBack)[0]
        usedfiles.append(filename[59:len(filename)]+', ')
        print('>---numpy.__version__ = {0}:'.format(np.__version__))
        print('>---We have within given range of coordinates:')
        print('                    {0} spectra in the East scan,'.format(np.shape(splitData.dataEast)[0]))
        print('                    {0} spectra in the Center scan,'.format(np.shape(splitData.dataCenter)[0]))
        print('                    {0} spectra in the West scan,'.format(np.shape(splitData.dataWest)[0]))
        print('                    {0} spectra in the back scan.\n'.format(np.shape(splitData.dataBack)[0]))
        print('These spectra correspond to the following coordinates:')
        print('+--------------------------+--------------------------+--------------------------+')
        print('|          East            |        Center            |          West            |')
        print('+--------------------------+--------------------------+--------------------------+')
        print('| Latitudes  |  Longitudes | Latitudes  |  Longitudes | Latitudes  |  Longitudes |')
        print('+--------------------------+--------------------------+--------------------------+')
        output = ''
        for i in range(len(latitudes)):
            if (not ma.array(latitudes[i]).all() is ma.masked):
                if (not ma.array(latitudes[i][0] is ma.masked)):
                    output += '| {0:9.5f}  |  {1:10.5f} '.format(latitudes[i,0],longitudes[i,0])
                else:
                    output += '|     --     |      --     '
                if (not ma.array(latitudes[i][1] is ma.masked)):
                    output += '| {0:9.5f}  |  {1:10.5f} '.format(latitudes[i,1],longitudes[i,1])
                else:
                    output += '|     --     |      --     '
                if (not ma.array(latitudes[i][2] is ma.masked)):
                    output += '| {0:9.5f}  |  {1:10.5f} |\n'.format(latitudes[i,2],longitudes[i,2])
                else:
                    output += '|     --     |      --     |\n'
        print(output[:-1])
        print('+--------------------------+--------------------------+--------------------------+')
        print('|       Backscan           |')
        print('+--------------------------+')
        print('| Latitudes  |  Longitudes |')
        print('+--------------------------+')
        for i in range(len(back_latitudes)):
            if (not ma.array(back_latitudes[i]).all() is ma.masked):
                if (not ma.array(back_latitudes[i][0] is ma.masked)):
                    print('| {0:9.5f}  |  {1:10.5f} |'.format(back_latitudes[i,0],back_longitudes[i,0]))
        print('+--------------------------+\n')
        if (do_irrad):
            idx += 1
        if (idx == 1):
            splitData.splitWavelength(inputData.wavelength,inputData.spectral_index[0],maskForward,inputData.startPixel,inputData.endPixel,inputData.irradiance_spectral_index)
            wavelengthIrradiance=splitData.wavelengthIrradiance
            del splitData.wavelengthIrradiance
        else:
            splitData.splitWavelength(inputData.wavelength,inputData.spectral_index[0],maskForward,inputData.startPixel,inputData.endPixel)
        splitData.splitWavelength(inputData.wavelength,inputData.back_spectral_index[0],maskBackward,inputData.back_startPixel,inputData.back_endPixel)
        endSplit = time.time()
        print('>---Time spent to split the data into East/Center/West/Backscan = {0:7.5f} seconds.---<'.format(endSplit-startSplit))
        
        # Interpolate data on the common grid
        startInterpol = time.time()
        print('>---Now we interpolate spectra onto a common wavelength grid.---<')
        if (do_irrad):
            splitData.interpolateData([0.0,0.0,0.0],do_irrad,wavelengthIrradiance)
        else:
            splitData.interpolateData(wavelength[0],do_irrad)
        endInterpol = time.time()
        print('>---Time spent to interpolate the data = {0:7.5f} seconds.---<'.format(endInterpol-startInterpol))
        
        # Calculate mean values of this input file and append the vectors to the master matrices
        startMean = time.time()
        print('>---Now we calculate the mean spectra for the current file.---<')
        splitData.calculateMean(splitData.interpolEast,splitData.interpolCenter,splitData.interpolWest,splitData.interpolBack)
        endMean = time.time()
        print('>---Time spent to average the data of the file = {0:7.5f} seconds.---<'.format(endMean-startMean))
        
        emptyEast = pa.DataFrame(splitData.meanEast).dropna().empty
        emptyCenter = pa.DataFrame(splitData.meanCenter).dropna().empty
        emptyWest = pa.DataFrame(splitData.meanWest).dropna().empty
        emptyBack = pa.DataFrame(splitData.meanBack).dropna().empty
        
        if (not emptyEast):
            masterEast.append(splitData.meanEast)
            nwavel=np.shape(splitData.meanEast)[0]
        if (not emptyCenter):
            masterCenter.append(splitData.meanCenter)
            nwavel=np.shape(splitData.meanCenter)[0]
        if (not emptyWest):
            masterWest.append(splitData.meanWest)
            nwavel=np.shape(splitData.meanWest)[0]
        if (not emptyBack):
            masterBack.append(splitData.meanBack)
            nwavel=np.shape(splitData.meanBack)[0]

        del splitData.interpolEast, splitData.interpolCenter, splitData.interpolWest, splitData.interpolBack
        del splitData.meanEast, splitData.meanCenter, splitData.meanWest, splitData.meanBack
        end = time.time()
        print('>---Total elapsed time for file = {0:7.5f} seconds.---<'.format(end-start))
        print('>---Continuing to the next file or exiting.---<\n')

    #sys.exit()
    print('+'+'='*(51)+'+')
    print('| Proceeding to calculate the total mean for today. |')
    print('+'+'='*(51)+'+\n')

    #print('netCDF.__version__ = ', netCDF4.__version__)
    emptyEast = pa.DataFrame(masterEast).T.dropna().empty
    emptyCenter = pa.DataFrame(masterCenter).T.dropna().empty
    emptyWest = pa.DataFrame(masterWest).T.dropna().empty
    emptyBack = pa.DataFrame(masterBack).T.dropna().empty

    if (emptyEast and emptyCenter and emptyWest and emptyBack):
        print('>---All of the master matrices (East/Center/West/backscan) are empty.---<')
        print('>---Exiting the program for today.---<')
        sys.exit()

    #Calculate total mean of the day from the master matrices
    startMean = time.time()
    print('>---Now we calculate the mean spectra for today.---<')
    splitData.calculateMean(masterEast,masterCenter,masterWest,masterBack)
    endMean = time.time()
    print('>---Time spent to average the data today = {0:7.5f} seconds.---<'.format(endMean-startMean))
    #np.savetxt(outfile,np.stack((wavelengthEval,splitData.meanEast,splitData.meanCenter,splitData.meanWest,splitData.meanBack)),delimiter='\n',newline='')

    if (emptyEast or emptyCenter or emptyWest or emptyBack):
        print('>---At least one of the master matrices (East/Center/West/backscan) is empty.---<')
        print('>---These will be filled with a _FillValue.---<')
    if emptyEast:
        splitData.meanEast=np.full(nwavel,fillValue)
    if emptyCenter:
        splitData.meanCenter=np.full(nwavel,fillValue)
    if emptyWest:
        splitData.meanWest=np.full(nwavel,fillValue)
    if emptyBack:
        splitData.meanBack=np.full(nwavel,fillValue)

    if (not do_irrad):
        wavelengthEval = np.linspace(wavelength[0][0],wavelength[0][1],(wavelength[0][1]-wavelength[0][0]) / wavelength[0][2]) 
        myframe=pa.DataFrame(np.stack((splitData.meanEast,splitData.meanCenter,splitData.meanWest,splitData.meanBack),axis=1),index=wavelengthEval,columns=['East','Center','West','Backscan']).fillna(0)
    else:
        myframe=pa.DataFrame(np.stack((splitData.meanEast,splitData.meanCenter,splitData.meanWest,splitData.meanBack),axis=1),index=wavelengthIrradiance,columns=['East','Center','West','Backscan']).fillna(0)
    if (do_apex):
        if (not do_irrad):
            write_apex(wavelengthEval,np.array(myframe.T),outputfile,numSpectra,latitude,longitude,usedfiles,fillValue)
        else:
            write_apex(wavelengthIrradiance,np.array(myframe.T),outputfile,numSpectra,latitude,longitude,usedfiles,fillValue)
    else:
        if (not do_irrad):
            write_text(wavelengthEval,np.array(myframe.T),outputfile)
        else:
            write_text(wavelengthIrradiance,np.array(myframe.T),outputfile)
    if do_plot:
        print('>---Plotting the results to the outputfile {0}.pdf.---<'.format(outputfile))
        if (not do_irrad):
            plot_result(wavelengthEval,splitData.meanEast,splitData.meanCenter,splitData.meanWest,splitData.meanBack,outputfile)
        else:
            plot_result(wavelengthIrradiance,splitData.meanEast,splitData.meanCenter,splitData.meanWest,splitData.meanBack,outputfile)

    #print(myframe.loc[myframe.all(axis=1)])
    #pa.DataFrame.to_dense(pa.DataFrame(myframe).isna().any(axis=1))
    #print(np.stack((wavelengthEval,splitData.meanEast,splitData.meanCenter,splitData.meanWest,splitData.meanBack)).shape)
    #print(np.stack((wavelengthEval,splitData.meanEast,splitData.meanCenter,splitData.meanWest,splitData.meanBack),axis=1).shape)
