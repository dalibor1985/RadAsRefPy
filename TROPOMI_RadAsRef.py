#!/bira-iasb/softs/18/py36/bin/python3.6

import os
import sys
import time
import numpy as np
import numpy.ma as ma
from netCDF4 import Dataset
#import netCDF4
#from tabulate import tabulate
from ReadNC import TROPOMIData
from DataAnalysis import AnalyseTROPOMI

def write_netcdf4(wavelengths,results,outfile,band,lat,lon,inputfiles,fillValue):
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
    dataset = Dataset(outfile, 'w',  format='NETCDF4')
    dataset.title = 'Reference file created by RAD_AS_REF for QDOAS Python script.'
    dataset.description = 'Radiance as reference file for QDOAS based on average radiances.'
    dataset.created = time.ctime(time.time())
    dataset.lat_bound = lat[0],lat[1]
    dataset.lon_bound = longitude[0],longitude[1]
    dataset.inputfiles = inputfiles
    dataset.measurement_date = "{0}/{1}/{2}".format(day,month,year)
    dataset.fillValue = fillValue
    obsName = band+'_RADIANCE/STANDARD_MODE/OBSERVATIONS'
    insName = band+'_RADIANCE/STANDARD_MODE/INSTRUMENT'
    observations = dataset.createGroup(obsName)
    instrument = dataset.createGroup(insName)
#    pixel = dataset.createDimension('pixel', np.shape(results)[0])
#    spectral_channel = dataset.createDimension('spectral_channel', np.shape(results)[1])
    reference_radiance = observations.createVariable('reference_radiance', np.float64, ('pixel','spectral_channel'))
    reference_wavelength = instrument.createVariable('reference_wavelength', np.float64, ('pixel','spectral_channel'))
    use_row = dataset.createVariable('use_row', np.int8, ('pixel',))
    reference_radiance[:] = results
    reference_wavelength[:] = wavelengths
    use_row[:] = np.where(np.all(np.isfinite(reference_radiance),axis=1),1,0)
    dataset.close()

def write_apex(wavelengths,results,outfile,lat,lon,inputfiles,fillValue):
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
    dataset.lat_bound = lat[0],lat[1]
    dataset.lon_bound = longitude[0],longitude[1]
    dataset.inputfiles = inputfiles
    dataset.measurement_date = "{0}/{1}/{2}".format(day,month,year)
    dataset.fillValue = fillValue
#    col_dim = dataset.createDimension('col_dim', np.shape(results)[0])
#    spectral_dim = dataset.createDimension('spectral_dim', np.shape(results)[1])
    reference_radiance = dataset.createVariable('reference_radiance', np.float64, ('col_dim','spectral_dim'))
    reference_wavelength = dataset.createVariable('reference_wavelength', np.float64, ('col_dim','spectral_dim'))
    use_row = dataset.createVariable('use_row', np.int8, ('col_dim',))
    reference_radiance[:] = results
    reference_wavelength[:] = wavelengths
    use_row[:] = np.where(np.all(np.isfinite(reference_radiance),axis=1),1,0)
    dataset.close()

def transformLongitudes(inputLongitudes,inputMask):
    assert(np.shape(inputLongitudes) == np.shape(inputMask))
    longitudes = ma.array(inputLongitudes,mask=inputMask)
    for i in range(len(longitudes)):
        for j in range(len(longitudes[i])):
            if (longitudes[i][j] < 0.0 and longitudes[i][j] > -180.):
                longitudes[i][j] = longitudes[i][j] + 360.0
    return longitudes

def mainTROPOMI(inputfiles,outputfile,latitude,longitude,longs,do_apex,*wavelength):
    # start the main loop over NetCDF4 files
    print('\n'+'%'*(48))
    print('%% RadAsRef.py script TROPOMI mode for QDOAS. %%')
    print('%'*(48))
    np.set_printoptions(threshold=np.inf,precision=5)
    np.warnings.filterwarnings('ignore')

    masterDay=[]
    idx = 0
    numSpectra=np.full(450,0)
    usedfiles = []
    for File in range(len(inputfiles)):
        startFile = time.time()
        filename=inputfiles[File]
        if (not filename.endswith('.nc')):
            filename += '.nc'
        if (not os.path.isfile(filename)):
            print('\n+'+'='*(30 + len(filename))+'+')
            print('| ERROR: file {0} does not exist! |'.format(filename))
            print('+'+'='*(30 + len(filename))+'+\n')
            sys.exit()
        else:
            print('\n+'+'='*(23 + len(filename))+'+')
            print('| Opening NetCDF4 file {0} |'.format(filename))
            print('+'+'='*(23 + len(filename))+'+\n')
        startRead = time.time()
        inputData = TROPOMIData()
        inputData.readCoordinates(filename)
        endRead = time.time()
        print('>---Time spent to read the coordinates = {0:7.5f} seconds.---<'.format(endRead-startRead))
        averageData=AnalyseTROPOMI()
        startMask = time.time()
        averageData.checkCoordinates(inputData.latitude,inputData.longitude,latitude,longitude,longs)
        endMask = time.time()
        print('>---Time spent to create the mask = {0:7.5f} seconds.---<'.format(endMask-startMask))
        if (np.all(averageData.mask) == True):
            print('''All pixels found outside of the requested coordinates.
                     Proceeding to the following file or leaving the loop.''')
            continue
        inputData = TROPOMIData()
        startRead = time.time()
        inputData.readRadiances(filename)
        endRead = time.time()
        print('>---Time spent to read the radiances = {0:7.5f} seconds.---<'.format(endRead-startRead))
        if (idx == 0):
            startRead = time.time()
            wavelengths=inputData.readWavelengths(filename)
            endRead = time.time()
            idx += 1
            print('>---Time spent to read the wavelengths = {0:7.5f} seconds.---<'.format(endRead-startRead))
        startChoose = time.time()
        averageData.chooseRadiances(inputData.radiance)
        fillValue=inputData._FillValue
        endChoose = time.time()
        print('>---Time spent to choose pixels = {0:7.5f} seconds.---<'.format(endChoose-startChoose))
        usedfiles.append(os.path.split(filename)[1])
        masterDay.append(np.delete(averageData.radiance,np.where(np.all(np.isnan(averageData.radiance),axis=(1,2))),axis=0))
        #continue
        #startAverage = time.time()
        #averageData.calculateMean()
        #endAverage = time.time()
        #print('>---Time spent to calculate the mean = {0:7.5f} seconds.---<'.format(endAverage-startAverage))
        #masterDay.append(averageData.meanFile[0])
        endFile = time.time()
        print('>---Time spent calculate the mean of the NetCDF4 file = {0:7.5f} seconds.---<'.format(endFile-startFile))

    print("len(masterDay) = ",len(masterDay))
    print("masterDay[0].shape = ",masterDay[0].shape)
    print("masterDay[1].shape = ",masterDay[1].shape)
    print("masterDay[2].shape = ",masterDay[2].shape)
    #try:
    print("np.concatenate((masterDay[0],masterDay[1]),axis=0): ",np.concatenate((masterDay[0],masterDay[1]),axis=0).shape)
    print("np.vstack(masterDay): ",np.vstack((masterDay))).shape
    #except MemoryError:
    #    print("MEMORY ERROR: lack of memory")
    #    return
    print("Exiting TROPOMI")
    return
    print('\n+'+'='*(51)+'+')
    print('| Proceeding to calculate the total mean for today. |')
    print('+'+'='*(51)+'+\n')

    try:
        startAverage = time.time()
        averageData.calculateMeanDay(np.dstack(masterDay))
        endAverage = time.time()
        print('>---Time spent to calculate the daily mean = {0:7.5f} seconds.---<'.format(endAverage-startAverage))
    except ValueError:
        print('No coordinates for this date fall within the given range.')
        return
    
    if (len(wavelength) > 0):
        print('\n+'+'='*(53)+'+')
        print('| Proceeding to interpolate the total mean for today. |')
        print('+'+'='*(53)+'+\n')
        startInterpolate = time.time()
        averageData.interpolateMean(wavelengths,wavelength[0][2])
        endInterpolate = time.time()
        print('>---Time spent to interpolate the daily mean = {0:7.5f} seconds.---<'.format(endInterpolate-startInterpolate))
    del masterDay
    if (do_apex):
        if (len(wavelength) > 0):
            write_apex(averageData.wavelengthEval,averageData.meanInterpolated,outputfile,latitude,longitude,usedfiles,fillValue)
        elif (len(wavelength) == 0):
            write_apex(wavelengths,averageData.meanDay,outputfile,latitude,longitude,usedfiles,fillValue)
    else: 
        if ('BD3' in inputfiles[0]):
            band = 'BAND3'
        if ('BD4' in inputfiles[0]):
            band = 'BAND4'
        if (len(wavelength) > 0):
            write_netcdf4(averageData.wavelengthEval,averageData.meanInterpolated,outputfile,band,latitude,longitude,usedfiles,fillValue)
        elif (len(wavelength) == 0):
            write_netcdf4(wavelengths,averageData.meanDay,outputfile,band,latitude,longitude,usedfiles,fillValue)

