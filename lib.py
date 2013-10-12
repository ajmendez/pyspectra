# Library functions
import os
import copy

import numpy as np
import pylab
import pyfits
from scipy import signal

import pysurvey.math
import pysurvey.plot
import pysurvey.file

print("Loading Library functions")
# http://skyserver.sdss3.org/dr9/en/tools/quicklook/quickobj.asp?plate=0282&mjd=51636&fiber=01


def lazyprop(fn):
    attr_name = '_lazy_' + fn.__name__
    @property
    def _lazyprop(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)
    return _lazyprop


def Show(fn):
    def _show(self, *args, **kwargs):
        noplot = kwargs.pop('noplot', False)
        quiet = kwargs.pop('quiet', False)
        tmp = fn(self, *args, **kwargs)
        if not quiet:
            self.info()
        if not noplot:
            self.plot(**kwargs)
        return tmp
    return _show

def Copy(fn):
    def _copy(self, *args, **kwargs):
        out = copy.copy(self)
        return fn(out)
    return _copy
        


class Spectra(object):
    '''Spectra is a 1d data structure'''
    @Show
    def __init__(self, filename, ext=0, index=0):
        '''Load the file -- for now only a fits file.'''
        self.filename = filename
        self.ext = ext
        self.index = index
        self.simplefilename = os.path.basename(self.filename)
        self.name = self.simplefilename
        self.action = 'Initialized'
        self.loaded = False
        self.loadFits()
    
    def __repr__(self):
        '''A simple representation of an object is to just give it a name'''
        return 'Spectra Object for {}'.format(self.name)
    
    @lazyprop
    def _wunit(self):
        '''Wavelength units'''
        return self.header['WAT1_001'].split('=')[-1]
    
    @lazyprop
    def _wrange(self):
        ii = np.where(self.flux > 0)
        xr = [np.min(self.wave[ii]), np.max(self.wave[ii])]
        return pysurvey.math.embiggen(xr,p=0.05)
    
    @lazyprop
    def _funit(self):
        return self.header['BUNIT']
    
    @lazyprop
    def _frange(self):
        yr = [0, np.max(self.flux)]
        return pysurvey.math.embiggen(yr,p=0.10, mode='upper')
        
    
    def info(self):
        '''Read some information from the header to inform the user'''
        print '''{}: {}  [{}]
 | Wavelength Coverage: [{:0.2f}, {:0.2f}] {}
 | Flux Units: {}'''.format(
     self.action, self.name, 
     # self.header['RA'], self.header['DEC'], 
     self.simplefilename,
     np.min(self.wave), np.max(self.wave), self._wunit,
     self._funit,
     )
    
    def loadFits(self):
        '''Grab the data'''
        self.header = pyfits.getheader(self.filename)
        self.data = pyfits.getdata(self.filename, self.ext)
        self.name = self.header['NAME']
        self.action = 'Loaded'
        
        self.wave = 10**(self.header['COEFF0']+
                         self.header['COEFF1']*
                         np.arange(0, self.header['NAXIS1']))
        self.flux = self.data[self.index]
        self.raw = copy.copy(self.flux)
    
    def plot(self, **kwargs):
        '''make the plot of the data'''
        tmp = {
            'xlabel':'Wavelength [{}]'.format(self._wunit),
            'xr': self._wrange,
            'ylabel': 'Flux [{}]'.format(self._funit),
            # 'yr': self.frange(),
        }
        tmp.update(kwargs)
        
        pylab.figure()
        pysurvey.plot.setup(**tmp)
        pylab.plot(self.wave, self.flux)
    
    @Copy
    @Show
    def __sub__(self,other):
        '''Subtraction of two Spectra Objects
        TODO: add @show @copy wrappers to simplify'''
        # out = copy.copy(self)
        out = self
        out.action = 'Subtracted'
        out.flux = out.flux - other.flux
        # out.info()
        # out.plot()
        return out
    
    @Copy
    @Show
    def getSmooth(self, smoothlen, name='boxcar'):
        '''Returns a smoothed spectra given a smoothing length (smoothlen).
        mode = [boxcar], triang, blackman, hamming, hann, bartlett, flattop, parzen, bohman, blackmanharris, nuttall, barthann, ...
        see scipy.signal.get_window()
        
        '''
        # out = copy.copy(self)
        out = self
        window = signal.get_window(name,smoothlen)
        out.flux = signal.convolve(self.flux, window/np.sum(window), mode='same')
        out.action = 'Smoothed'
        # out.info()
        # out.plot()
        return out
    
    
    def getLine(self, linename='Lick_Mg2'):
        '''Get a fit line from the fits file'''
        linetable = pysurvey.file.Cat(self.filename, 5)
        ii = np.where(linetable['name'] == linename)
        line = linetable[ii]
        r = [line.waveMin,
             line.waveMax]
        tmp = {
            'xr': pysurvey.math.embiggen(r, p=2)
        }
        print ''.join(line.__str__().splitlines(True)[1:])
        self.plot( **tmp)
        pysurvey.plot.line(x=r)
    
    
    

def loadSpectra(filename):
    '''Load spectra from a given filename'''
    if '.fits' not in filename:
        raise ValueError('File: {} should be a .fits/.fits.gz file'.format(filename))
    return Spectra(filename)
    