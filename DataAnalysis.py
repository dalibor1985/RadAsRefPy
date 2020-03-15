#!/bira-iasb/softs/18/py36/bin/python3.6

import sys
import numpy as np
import numpy.ma as ma
import pandas as pa
import socket
from scipy import interpolate

class AnalyseGOME():

    def __init__(self):
        self.dataEast=[]
        self.dataCenter=[]
        self.dataWest=[]
        self.dataBack=[]
        self.wavelengthEast=[]
        self.wavelengthCenter=[]
        self.wavelengthWest=[]
        self.wavelengthBack=[]
        self.interpolEast=[]
        self.interpolCenter=[]
        self.interpolWest=[]
        self.interpolBack=[]
        pass

    def splitSpectra(self,rawSpectra,inputMask):
        assert((np.shape(inputMask)[1] == 3) or (np.shape(inputMask)[1] == 1))
        print('socket.gethostname() = ', socket.gethostname())
        if (socket.gethostname() == "dalibornb"):
            rawSpectra[:,:,0] = ma.array(rawSpectra[:,:,0],mask=inputMask)
            for i in range(len(rawSpectra)):
                rawSpectra[i] = ma.array(ma.mask_rows(rawSpectra[i],axis=1))
        if (np.shape(inputMask)[1] == 3): # masks for forward
            if (socket.gethostname() == "dalibornb"):
                self.dataEast = rawSpectra[:,0,:][~np.all(np.isnan(rawSpectra[:,0,:]),axis=1)]
                self.dataCenter = rawSpectra[:,1,:][~np.all(np.isnan(rawSpectra[:,1,:]),axis=1)]
                self.dataWest = rawSpectra[:,2,:][~np.all(np.isnan(rawSpectra[:,2,:]),axis=1)]
            else:
                for i in range(len(rawSpectra)):
                    if (inputMask[i][0].all() == False):
                        self.dataEast.append(rawSpectra[i,0,:])
                    if (inputMask[i][1].all() == False):
                        self.dataCenter.append(rawSpectra[i,1,:])
                    if (inputMask[i][2].all() == False):
                        self.dataWest.append(rawSpectra[i,2,:])

        elif (np.shape(inputMask)[1] == 1): # for backscan
            if (socket.gethostname() == "dalibornb"):
                self.dataBack = rawSpectra[:,0,:][~np.all(np.isnan(rawSpectra[:,0,:]),axis=1)]
            else:
                for i in range(len(rawSpectra)):
                    if (inputMask[i].all() == False):
                        self.dataBack.append(rawSpectra[i,0,:])

        return

    def splitWavelength(self,rawWavelength,spectralIndex,inputMask,startPix,endPix,*irradianceSpectralIndex):
        assert((np.shape(spectralIndex)[1] == 3) | (np.shape(spectralIndex)[1] == 1))
        assert(np.shape(spectralIndex) == np.shape(inputMask))
        if (np.shape(spectralIndex)[1] == 3): # mask for forward
            for i in range(np.shape(spectralIndex)[0]):
                if (spectralIndex[i][0] >= 0 and spectralIndex[i][0] < 8 and inputMask[i][0] == False):
                    self.wavelengthEast.append(rawWavelength[spectralIndex[i][0]][1][startPix:endPix+1])
                if ( spectralIndex[i][1] >= 0 and spectralIndex[i][1] < 8 and inputMask[i][1] == False):
                    self.wavelengthCenter.append(rawWavelength[spectralIndex[i][1]][1][startPix:endPix+1])
                if ( spectralIndex[i][2] >= 0 and spectralIndex[i][2] < 8 and inputMask[i][2] == False):
                    self.wavelengthWest.append(rawWavelength[spectralIndex[i][2]][1][startPix:endPix+1])
        elif (np.shape(spectralIndex)[1] == 1): # for backscan
            for i in range(np.shape(spectralIndex)[0]):
                if (spectralIndex[i][0] >= 0 and spectralIndex[i][0] < 8 and inputMask[i][0] == False):
                    self.wavelengthBack.append(rawWavelength[spectralIndex[i][0]][1][startPix:endPix+1])
        if (len(irradianceSpectralIndex) > 0):
            self.wavelengthIrradiance=rawWavelength[irradianceSpectralIndex][1][:]
        return

#    def createMask(self):
#        self.maskData=(pa.DataFrame(self.dataEast).isna().any(axis=1) | pa.DataFrame(self.dataCenter).isna().any(axis=1) | pa.DataFrame(self.dataWest).isna().any(axis=1))

    def maskGeodata(self,latitudes,longitudes,lat_bounds,lon_bounds,longs):
        assert(np.shape(latitudes) == np.shape(longitudes))
        np.warnings.filterwarnings('ignore')
        maskLatitudes = ((latitudes<lat_bounds[0]) | (latitudes>lat_bounds[1])) | np.isnan(latitudes)
        if (longs == 1):
            maskLongitudes = ((longitudes<lon_bounds[0]) & (longitudes>lon_bounds[1])) | np.isnan(longitudes)
        if (longs == 2):
            maskLongitudes = ((longitudes>lon_bounds[0]) & (longitudes<lon_bounds[1])) | np.isnan(longitudes)
        return maskLatitudes | maskLongitudes

    def interpolateData(self,wavelength,do_irrad,*wavelengthIrradiance):
        emptyEast = pa.DataFrame(self.dataEast).dropna().empty
        emptyCenter = pa.DataFrame(self.dataCenter).dropna().empty
        emptyWest = pa.DataFrame(self.dataWest).dropna().empty
        emptyBack = pa.DataFrame(self.dataBack).dropna().empty
 
        if (emptyEast and emptyCenter and emptyWest and emptyBack):
            print('all is empty')
            return 
        if (do_irrad):
            wavelengthEval=wavelengthIrradiance
        else:
            npoints = (wavelength[1]-wavelength[0]) / wavelength[2]
            wavelengthEval = np.linspace(wavelength[0],wavelength[1],npoints) 
        
        if (not emptyEast):
            for i in range(len(self.wavelengthEast)):
                self.interpolEast.append(interpolate.griddata(self.wavelengthEast[i],self.dataEast[i],wavelengthEval,method='cubic',rescale=False))
        if (not emptyCenter):
            for i in range(len(self.wavelengthCenter)):
                self.interpolCenter.append(interpolate.griddata(self.wavelengthCenter[i],self.dataCenter[i],wavelengthEval,method='cubic',rescale=False))
        if (not emptyWest):
            for i in range(len(self.wavelengthWest)):
                self.interpolWest.append(interpolate.griddata(self.wavelengthWest[i],self.dataWest[i],wavelengthEval,method='cubic',rescale=False))
        if (not emptyBack):
            for i in range(len(self.wavelengthBack)):
                self.interpolBack.append(interpolate.griddata(self.wavelengthBack[i],self.dataBack[i],wavelengthEval,method='cubic',rescale=False))
        return

    def calculateMean(self,dataEast,dataCenter,dataWest,dataBack):
        if (not pa.DataFrame(dataEast).T.dropna().empty):
            self.meanEast=np.ma.average(np.ma.array(dataEast),axis=0)
        else:
            self.meanEast=[]
        if (not pa.DataFrame(dataCenter).T.dropna().empty):
            self.meanCenter=np.ma.average(np.ma.array(dataCenter),axis=0)
        else:
            self.meanCenter=[]
        if (not pa.DataFrame(dataWest).T.dropna().empty):
            self.meanWest=np.ma.average(np.ma.array(dataWest),axis=0)
        else:
            self.meanWest=[]
        if (not pa.DataFrame(dataBack).T.dropna().empty):
            self.meanBack=np.ma.average(np.ma.array(dataBack),axis=0)
        else:
            self.meanBack=[]
        return

class AnalyseTROPOMI():

    def __init__(self):
        pass

    def checkCoordinates(self,latitudes,longitudes,lat_bounds,lon_bounds,longs):
        assert(np.shape(latitudes) == np.shape(longitudes))
        maskLatitudes = (latitudes<lat_bounds[0]) | (latitudes>lat_bounds[1])
        if (longs == 1):
            maskLongitudes = (longitudes<lon_bounds[0]) & (longitudes>lon_bounds[1])
        if (longs == 2):
            maskLongitudes = (longitudes>lon_bounds[0]) & (longitudes<lon_bounds[1])
        self.mask=(maskLatitudes | maskLongitudes)

    def chooseRadiances(self,rawSpectra):
        assert(np.shape(rawSpectra)[:2] == np.shape(self.mask))
        compoundMask=np.dstack([self.mask]*np.shape(rawSpectra)[2])
        #try:
        self.radiance=np.where(compoundMask==True,np.nan,rawSpectra)
        #except MemoryError:
        #    print("There is a memory error!")
        #    sys.exit()

    def calculateMean(self):
        print("self.radiance.shape = ",self.radiance.shape)
        self.meanFile=np.nanmean(self.radiance,axis=0,dtype='float64',keepdims=True)
        print("self.meanFile.shape = ",self.meanFile.shape)

    def calculateMeanDay(self,masterDay):
        self.meanDay=np.nanmean(masterDay,axis=2,dtype='float64')

    def interpolateMean(self,wavelengths,step):
        npoints=int(np.amax(np.ceil((wavelengths[:,-1]-wavelengths[:,0])/step)))
        print("npoints = ",npoints)
        #try:
        self.meanInterpolated=np.empty((np.shape(self.meanDay)[0],npoints))
        self.wavelengthEval=np.empty((np.shape(self.meanDay)[0],npoints))
        for i in range(np.shape(self.meanDay)[0]):
            self.wavelengthEval[i] = np.linspace(wavelengths[i,0],wavelengths[i,-1],npoints) 
            self.meanInterpolated[i]=interpolate.griddata(wavelengths[i],self.meanDay[i],self.wavelengthEval[i],method='cubic',rescale=False)
        #except MemoryError:
        #    print("There is a memory error!")
        #    sys.exit()
