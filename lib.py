# Library functions
import os
import copy

import numpy as np
import pylab
import pyfits
from scipy import signal

import pysurvey.math
import pysurvey.plot

print("Loading Library functions")
# http://skyserver.sdss3.org/dr9/en/tools/quicklook/quickobj.asp?plate=0282&mjd=51636&fiber=01


class Spectra(object):
    def __init__(self, filename):
        self.filename = filename
        self.simplefilename = os.path.basename(self.filename)
        self.name = self.simplefilename
        self.action = 'Initialized'
        self.loaded = False
        self.loadFits()
        self.info()
        self.plot()
    
    def __repr__(self):
        return 'Spectra Object for {}'.format(self.name)
    
    def name(self):
        '''Return a simple name for the object'''
        return self.name
    
    def wunit(self):
        '''Wavelength units'''
        return self.header['WAT1_001'].split('=')[-1]
    def wrange(self):
        ii = np.where(self.flux > 0)
        xr = [np.min(self.wave[ii]), np.max(self.wave[ii])]
        return pysurvey.math.embiggen(xr,p=0.05)
        
    def funit(self):
        return self.header['BUNIT']
    def frange(self):
        yr = [0, np.max(self.flux)]
        return pysurvey.math.embiggen(yr,p=0.10, mode='upper')
        
    
    def info(self):
        '''Read some information from the header to inform the user'''
        print '''{}: {} [{}]
 | RA/DEC: {}, {}
 | Wavelength Coverage: [{:0.2f}, {:0.2f}] {}
 | Flux Units: {}'''.format(
     self.action, self.name, self.simplefilename,
     self.header['RA'], self.header['DEC'],
     np.min(self.wave), np.max(self.wave), self.wunit(),
     self.funit(),
     )
    
    def loadFits(self):
        '''Grab the data'''
        self.header = pyfits.getheader(self.filename)
        self.data = pyfits.getdata(self.filename)
        self.name = self.header['NAME']
        self.action = 'Loaded'
        
        self.wave = 10**(self.header['COEFF0']+
                         self.header['COEFF1']*
                         np.arange(0, self.header['NAXIS1']))
        self.flux = self.data[0]
    
    def plot(self, **kwargs):
        
        tmp = {
            'xlabel':'Wavelength [{}]'.format(self.wunit()),
            'xr': self.wrange(),
            'ylabel': 'Flux [{}]'.format(self.funit()),
            # 'yr': self.frange(),
        }
        tmp.update(kwargs)
        
        pylab.figure()
        pysurvey.plot.setup(**tmp)
        pylab.plot(self.wave, self.flux)
    
    
    def __sub__(self,other):
        '''Subtraction of two Spectra Objects
        TODO: add @show @copy wrappers to simplify'''
        out = copy.copy(self)
        out.action = 'Subtracted'
        out.flux = out.flux - other.flux
        out.info()
        out.plot()
        return out
    
    def getSmooth(self, smoothlen, name='boxcar'):
        '''Returns a smoothed spectra given a smoothing length (smoothlen).
        mode = [boxcar], triang, blackman, hamming, hann, bartlett, flattop, parzen, bohman, blackmanharris, nuttall, barthann, ...
        see scipy.signal.get_window()
        
        '''
        out = copy.copy(self)
        window = signal.get_window(name,smoothlen)
        out.flux = signal.convolve(self.flux, window/np.sum(window), mode='same')
        out.action = 'Smoothed'
        out.info()
        out.plot()
        return out
    
    def fitEV(self, wcenter, wrange):
        pass
        

    
    
    

def loadSpectra(filename):
    '''Load spectra from a given filename'''
    if '.fits' not in filename:
        raise ValueError('File: {} should be a .fits/.fits.gz file'.format(filename))
    return Spectra(filename)
    