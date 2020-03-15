#!/bira-iasb/softs/18/py36/bin/python3.6

import os
import sys
import math
import datetime as dt
import numpy as np
import numpy.ma as ma
import pandas as pa
import string
import argparse
from netCDF4 import Dataset

class GOMEData():

    def __init__(self):
        pass

    def read(self,filename):
        self.dataset = Dataset(filename,'r')
        self.empty_nadir = False
        try:
            self.calibration = self.dataset.groups['CALIBRATION']
            self.irradiance = self.dataset.groups['IRRADIANCE']
            self.observations = self.dataset.groups['MODE_NADIR'].groups['BAND_2B'].groups['OBSERVATIONS']
            self.geoData = self.dataset.groups['MODE_NADIR'].groups['BAND_2B'].groups['GEODATA']
            self.back_observations = self.dataset.groups['MODE_NADIR_BACKSCAN'].groups['BAND_2B'].groups['OBSERVATIONS']
            self.back_geoData = self.dataset.groups['MODE_NADIR_BACKSCAN'].groups['BAND_2B'].groups['GEODATA']
        except KeyError:
            self.empty_nadir = True
            return
        self.startPixel = self.dataset.groups['MODE_NADIR'].groups['BAND_2B'].variables['start_pixel'][:]
        self.endPixel = self.dataset.groups['MODE_NADIR'].groups['BAND_2B'].variables['end_pixel'][:]
        self.back_startPixel = self.dataset.groups['MODE_NADIR_BACKSCAN'].groups['BAND_2B'].variables['start_pixel'][:]
        self.back_endPixel = self.dataset.groups['MODE_NADIR_BACKSCAN'].groups['BAND_2B'].variables['end_pixel'][:]
        self.FillValue = self.observations.variables['radiance']._FillValue
        self.getTemperature()
        self.getWavelength()
        self.getSpectralIndex()
        self.getRadiance()
        self.getGeodata()
        self.dataset.close()

    def getTemperature(self):
        self.temperature=GOMEData.makeMasked(self.calibration.variables['temperature'][:])

    def getWavelength(self):
        self.wavelength=GOMEData.makeMasked(self.calibration.variables['wavelength'][:])

    def getSpectralIndex(self):
        self.spectral_index = GOMEData.makeMasked(self.observations.variables['spectral_index'][:])
        self.back_spectral_index = GOMEData.makeMasked(self.back_observations.variables['spectral_index'][:])
        self.spectral_index = self.spectral_index.astype(int)
        self.back_spectral_index = self.back_spectral_index.astype(int)
        self.irradiance_spectral_index=self.irradiance.variables["spectral_index"][:]

    def getRadiance(self):
        self.radiance = GOMEData.makeMasked(self.observations.variables['radiance'][:])
        self.back_radiance = GOMEData.makeMasked(self.back_observations.variables['radiance'][:])

    def getGeodata(self):
        self.latitudes = GOMEData.makeMasked(self.geoData.variables['latitude'][:][:])
        self.longitudes = GOMEData.makeMasked(self.geoData.variables['longitude'][:][:])
        self.back_latitudes = GOMEData.makeMasked(self.back_geoData.variables['latitude'][:][:])
        self.back_longitudes = GOMEData.makeMasked(self.back_geoData.variables['longitude'][:][:])

    @staticmethod
    def makeMasked(arr):
        ''' when using netcdf4 from python it will read data without masked values in it as numpy array, so far so good. When reading netcdf4 with masked values, we see that netcdf4 lib create a maskedarray, we do not want this, so in this function we change the masked values in a mask array to numpy arrays with nan, with makes our life a lot easier. The only reason why WE would use masked arrays is when we want to remove from time to time the mask, but this is not really the case.'''
        if type(arr) == ma.core.MaskedArray:
            return np.where(arr.mask==True,np.nan,arr )
        else: return arr

class TROPOMIData():

    def __init__(self):
        pass

    def readCoordinates(self,filename):
        self.dataset = Dataset(filename,'r') 
        try:
            if ('BD3' in filename):
                geoData=self.dataset.groups['BAND3_RADIANCE'].groups['STANDARD_MODE'].groups['GEODATA']
            elif ('BD4' in filename):
                geoData=self.dataset.groups['BAND4_RADIANCE'].groups['STANDARD_MODE'].groups['GEODATA']
            else:
                print('For now we are only reading bands 3 and 4!')
                print('Change your input!')
                sys.exit()
        except KeyError:
            print('Group GEODATA not found!')
            return
        try:
            self.latitude=TROPOMIData.makeMasked(geoData.variables['latitude'][0])
            self.longitude=TROPOMIData.makeMasked(geoData.variables['longitude'][0])
        except KeyError:
            print('Variable latitute and/or longitude not found!')
            return
        self.dataset.close()

    def readRadiances(self,filename):
        self.dataset = Dataset(filename,'r') 
        try:
            if ('BD3' in filename):
                observations=self.dataset.groups['BAND3_RADIANCE'].groups['STANDARD_MODE'].groups['OBSERVATIONS']
            elif ('BD4' in filename):
                observations=self.dataset.groups['BAND4_RADIANCE'].groups['STANDARD_MODE'].groups['OBSERVATIONS']
            else:
                print('For now we are only reading bands 3 and 4!')
                sys.exit()
        except KeyError:
            print('Group OBSERVATIONS not found!')
            return
        try:
            self.radiance=TROPOMIData.makeMasked(observations.variables['radiance'][0])
            self._FillValue=observations.variables['radiance']._FillValue
        except KeyError:
            print('Variable radiance not found!')
            return
        #except MemoryError:
        #    print('There is a memory error!')
        #    return
        self.dataset.close()

    def readWavelengths(self,filename):
        self.dataset = Dataset(filename,'r') 
        try:
            if ('BD3' in filename):
                instrument=self.dataset.groups['BAND3_RADIANCE'].groups['STANDARD_MODE'].groups['INSTRUMENT']
            elif ('BD4' in filename):
                instrument=self.dataset.groups['BAND4_RADIANCE'].groups['STANDARD_MODE'].groups['INSTRUMENT']
            else:
                print('For now we are only reading bands 3 and 4!')
                sys.exit()
        except KeyError:
            print('Group INSTRUMENT not found!')
            return
        try:
            wavelength=TROPOMIData.makeMasked(instrument.variables['nominal_wavelength'][0])
        except KeyError:
            print('Variable nominal_wavelength not found!')
            return
        self.dataset.close()
        return wavelength

    @staticmethod
    def makeMasked(arr):
        ''' when using netcdf4 from python it will read data without masked values in it as numpy array, so far so good. When reading netcdf4 with masked values, we see that netcdf4 lib create a maskedarray, we do not want this, so in this function we change the masked values in a mask array to numpy arrays with nan, with makes our life a lot easier. The only reason why WE would use masked arrays is when we want to remove from time to time the mask, but this is not really the case.'''
        if type(arr) == ma.core.MaskedArray:
            return np.where(arr.mask==True,np.nan,arr)
        else: return arr

