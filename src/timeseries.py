import os
import numpy as np
import mtspec
import xarray as xr
from statsmodels.tsa.arima_process import ArmaProcess

from paths import path_ocn_ctrl, path_ocn_rcp
from paths import path_atm_ctrl, path_atm_rcp
from paths import path_yrly_ctrl, path_yrly_rcp
from paths import rcpstr, spinup, CESM_filename
from paths import path_samoc, path_results
from filters import lowpass, chebychev, notch, deseasonalize  # old code imports those from here
from constants import abs_zero
from xr_integrate import xr_surf_mean, xr_zonal_mean
from xr_DataArrays import xr_AREA


class IterateOutputCESM:
    """ iterator over all CESM ctrl/rcp filenames
    automatically detects the last file
    
    example:
    >for year, month, filename in IterateOutputCESM('ocn', ctrl', 'monthly'):
    >    print(year, month, filename)
    """
    
    def __init__(self, domain, run, tavg, name=None):
        assert domain in ['ocn', 'ocn_rect', 'ocn_low', 'atm', 'ice']
        assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'pop']
        assert tavg in ['monthly', 'yrly']
        
        self.domain = domain
        self.run    = run
        self.tavg   = tavg
        self.stop   = False
        self.name   = name
        
        if tavg=='monthly':  self.month = 1
        elif tavg=='yrly':   self.month = 0
        
        if run=='ctrl':   self.year  =  100
        elif run=='rcp':  self.year  = 2000
        elif run=='pop':  self.year  =  125
        elif run=='lpd':  self.year  =  154
        elif run=='lpi':  
            if domain in ['ocn', 'ocn_low']:
                self.year  = 1600
            if domain=='atm':
                self.year  = 2876
            
    def file(self):
        if self.tavg=='monthly':
            filename = CESM_filename(self.domain, self.run, self.year, self.month)     
        elif self.tavg=='yrly':
            if self.domain=='ocn':
                if self.name==None:
                    raise ValueError('must provide (variables part of) name for yrly file')
                else:
                    filename = CESM_filename(self.domain, self.run, self.year, self.month, self.name)
            elif self.domain=='atm':
                if self.run in ['ctrl', 'rcp', 'lpi']:
                    filename = CESM_filename(self.domain, self.run, self.year, self.month, self.name)
                elif self.run in ['lpd']:  # yrly files are written out already
                    if self.name!=None:  print("name is ignored, as yearly files existed already")
                    filename = CESM_filename(domain=self.domain, run=self.run, y=self.year, m=self.month, name=None)
        if os.path.exists(filename)==False:
            self.stop = True
        return filename
    
    def __iter__(self):
        return self
    
    def __len__(self):
        length = 0
        year   = self.year
        month  = self.month
        while os.path.exists(self.file()):
            if self.tavg=='monthly':
                self.month += 1
                length +=1
                if self.month==13:
                    self.month = 1
                    self.year +=1
            elif self.tavg=='yrly':
                self.year +=1
                length +=1
        return length

    def __next__(self):
        new_file = self.file()
        y = self.year
        m = self.month
        if self.stop:
            raise StopIteration
        else:
            if self.tavg=='monthly':
                if self.month==12:
                    self.month = 1
                    self.year +=1
                else:
                    self.month += 1
                return y, m, new_file
            elif self.tavg=='yrly':
                self.year += 1
                return y, 0, new_file



def ncfile_list(domain, run, tavg, name=None):
    list_files = []
    for (y,m,s) in IterateOutputCESM(domain=domain, run=run, tavg=tavg, name=name):
        assert os.path.exists(s)
        list_files.append(s)
    return list_files



class TimeSeriesAnalysis(object):
    """ commonly used time series analysis functions 
    > spectrum
    > AR(1) fitting
    > MC spectral (uncertainty) estimates
    
    """
    
    def __init__(self, ts):
        """ time series analysis
        ts .. (np.ndarray) regularly sampled time series data
        """
        self.ts = ts
        assert type(self.ts)==np.ndarray
        self.len = len(self.ts)
    
    def spectrum(self, data=None, filter_type=None, filter_cutoff=None):
        """ multitaper spectrum """
        if data is None:  data = self.ts
        assert type(data)==np.ndarray
        if filter_type is not None:
            assert filter_type in ['lowpass', 'chebychev']
            assert type(filter_cutoff)==int
            assert filter_cutoff>1
            if filter_type=='lowpass':
                data = lowpass(data, filter_cutoff)
            elif filter_type=='chebychev':
                data = chebychev(data, filter_cutoff)
        
        spec, freq, jackknife, _, _ = mtspec.mtspec(
                data=data, delta=1., time_bandwidth=4,
                number_of_tapers=5, statistics=True)
        return (spec, freq, jackknife)
    
    def autocorrelation(self, n=1):
        """ calculates the first n lag autocorrelation coefficient """
        assert n<self.len
        n += 1  # zeroth element is lag-0 autocorrelation = 1
        acs = np.ones((n))
        for i in np.arange(1,n):
            acs[i] = np.corrcoef(self.ts[:-i],self.ts[i:])[0,1]
        return acs
        
    def mc_ar1(self, n=1000):
        """ Monte-Carlo AR(1) processes """
        phi = self.autocorrelation(n=1)[1]
        AR_object = ArmaProcess(np.array([1, -phi]), np.array([1]))
        mc = np.zeros((n, self.len))
        for i in range(n):
            mc[i,:] = AR_object.generate_sample(nsample=self.len)
            mc[i,:] *= np.std(self.ts)/np.std(mc[i,:])
        return mc
        
    def mc_ar1_spectrum(self, n=1000, filter_type=None, filter_cutoff=None):
        """ calculates the MC avg spectrum and the 95% confidence interval """
        mc = self.mc_ar1()
        if filter_type is not None:
            assert filter_type in ['lowpass', 'chebychev']
            assert type(filter_cutoff)==int
            assert filter_cutoff>1
            if filter_type=='lowpass':
                mc = lowpass(mc.T, filter_cutoff).T
            elif filter_type=='chebychev':
                mc = chebychev(mc.T, filter_cutoff).T

        mc_spectra = np.zeros((n, int(self.len/2)+1))
        for i in range(n):
            (mc_spectra[i,:], freq, jk) = self.spectrum(data=mc[i,:])
        mc_spectrum = np.zeros((4, int(self.len/2)+1))
        mc_spectrum[0,:] = np.median(mc_spectra, axis=0)
        mc_spectrum[1,:] = freq
        mc_spectrum[2,:] = np.percentile(mc_spectra,  5, axis=0)
        mc_spectrum[3,:] = np.percentile(mc_spectra, 95, axis=0)
        return mc_spectrum


            

class TimeSeriesSynthesis(object):
    """ combines different time series """
    def __init__(self):
        self.ctrl = {'name':'ctrl'}
        self.rcp  = {'name':'rcp'}
        self.lpd  = {'name':'lpd'}
        self.lpi  = {'name':'lpi'}
        self.had  = {'name':'had'}
        
        self.runs = [self.ctrl, self.rcp, self.lpd, self.lpi, self.had]
        
        self.load_GMST()
        self.load_TOA()
        self.load_SST_indices()
        
        
    def load_TOA(self):
        return
        
    def load_GMST(self):
        for run in self.runs[:4]:
            name = run["name"]
            run['GMST'] = xr.open_dataset(f'{path_samoc}/GMST/GMST_yrly_{name}.nc', decode_times=False).GMST
        return
    
    def load_SST_indices(self):
        for run in self.runs:
            for index in ['AMO', 'SOM', 'TPI']:
                name = run["name"]
                run[index] = xr.open_dataarray(f'{path_samoc}/SST/{index}_{name}.nc', decode_times=False)
        return
        
    def plot_lead_lag(self):
        return
    
    def model_GMST_SST_indices(self):
        """ linear models of GMST as a function of SST indices """


