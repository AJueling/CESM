"""
contains 3 classes:
"""

import numpy as np
import mtspec
import xarray as xr
from statsmodels.tsa.arima_process import ArmaProcess

from filters import lowpass


class xrAnalysis(object):
    """ functions that work on 1D time series as well as on 3D fields
    > trends
    """
    
    def __init__(self):
        return
    
    def lintrend(self, x):
        """ linear trend timeseries of a timeseries """
        pf = np.polynomial.polynomial.polyfit(x.time, x, 1)
        lf = pf[1]*x.time + pf[0]
        return lf
    
    def quadtrend(self, x):
        """ quadratic trend timeseries of a timeseries """
        pf = np.polynomial.polynomial.polyfit(x.time, x, 2)
        lf = pf[2]*x.time**2 + pf[1]*x.time + pf[0]
        return lf
    
    def standard_deviation(self, x):
        """ calculated standard deviation
        can deal with NaNs in time dimension, e.g. in HadISST
        """    
        xstd = x.std(axis=0, skipna=True)
        xstd.name = 'standard_deviation'
        return xstd
    
    def autocorrelation(self, x):
        """ autocorrelation 
        can deal with NaNs in time dimension, e.g. in HadISST
        """    
        x, y = x[1:,:,:], x.shift(time=1)[1:,:,:]
        y = xr.where(np.isnan(x), np.nan, y)
        x = xr.where(np.isnan(y), np.nan, x)

        n     = xr.where(np.isnan(x), 0, 1).sum(axis=0)
        xmean = x.mean(axis=0, skipna=True)
        ymean = y.mean(axis=0, skipna=True)
        xstd  = x.std( axis=0, skipna=True)
        ystd  = y.std( axis=0, skipna=True)

        x -= xmean
        y -= ymean

        cov      = np.divide(np.nansum(np.multiply(x,y), axis=0), n)
        cor      = cov/(xstd*ystd)
        cor.name = 'autocorrelation'

        return cor
    
    def lag_linregress(self, x, y, dof_corr=1, lagx=0, lagy=0):
        """
        adapted from: https://hrishichandanpurkar.blogspot.com/2017/09/vectorized-functions-for-correlation.html
        Input: Two xr.Datarrays of any dimensions with the first dim being time. 
        Thus the input data could be a 1D time series, or for example, have three dimensions (time,lat,lon). 
        Datasets can be provied in any order, but note that the regression slope and intercept will be calculated
        for y with respect to x.
        Output: xr Dataset containing covariance, correlation, regression slope and intercept, p-value, and
        standard error on regression between the two datasets along their aligned time dimension.  
        Lag values can be assigned to either of the data, with lagx shifting x, and lagy shifting y, with the specified lag amount.
        dof_corr .. (0,1] correction factor for reduced degrees of freedom
        """ 
        assert dof_corr<=1 and dof_corr>0
        #1. Ensure that the data are properly aligned to each other.
        x,y = xr.align(x,y)

        #2. Add lag information if any, and shift the data accordingly
        if lagx!=0:
            #If x lags y by 1, x must be shifted 1 step backwards. 
            #But as the 'zero-th' value is nonexistant, xr assigns it as invalid (nan). Hence it needs to be dropped
            x   = x.shift(time = -lagx).dropna(dim='time')
            #Next important step is to re-align the two datasets so that y adjusts to the changed coordinates of x
            x,y = xr.align(x,y)

        if lagy!=0:
            y   = y.shift(time = -lagy).dropna(dim='time')
            x,y = xr.align(x,y)

        #3. Compute data length, mean and standard deviation along time axis for further use: 
        n     = x.shape[0]
        xmean = x.mean(axis=0)
        ymean = y.mean(axis=0)
        xstd  = x.std(axis=0)
        ystd  = y.std(axis=0)

        #4. Compute covariance along time axis
        cov   =  np.sum((x - xmean)*(y - ymean), axis=0)/(n)

        #5. Compute correlation along time axis
        cor   = cov/(xstd*ystd)

        #6. Compute regression slope and intercept:
        slope     = cov/(xstd**2)
        intercept = ymean - xmean*slope  

        #7. Compute P-value and standard error
        #Compute t-statistics
        tstats = cor*np.sqrt(n*dof_corr-2)/np.sqrt(1-cor**2)
        stderr = slope/tstats

        pval   = sp.stats.t.sf(tstats, n-2)  # *2 for t-tailed test
        pval   = xr.DataArray(pval, dims=cor.dims, coords=cor.coords)

        cov.name       = 'cov'
        cor.name       = 'cor'
        slope.name     = 'slope'
        intercept.name = 'intercept'
        pval.name      = 'pval'
        stderr.name    = 'stderr'

        ds = xr.merge([cov, cor, slope, intercept, pval, stderr])

        ds.attrs['first_year'] = int(y.time[0]/365)
        ds.attrs['last_year']  = int(y.time[-1]/365)
        ds.attrs['lagx'] = lagx
        ds.attrs['lagy'] = lagy

        return ds

    
class TimeSeriesAnalysis(xrAnalysis):
    """ functions to analyze single 1D xr time series  
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
    
    
    
class FieldAnalysis(xrAnalysis):
    """ functions to analyze single 3D xr fields
    > trend map
    > std map
    > autocorrelation map
    
    """
    def __init__(self, field):
        """
        field .. xr.DataArray
        """
        self.field = field
    
    def make_linear_trend_map(self, fn):
        """ """
        self.load_SST_data()
        return
        
    def make_standard_deviation_map(self, fn=None):
        ds = standard_deviation(self.field)
        if fn is not None:  ds.to_netcdf(fn)
        return ds
    
    def make_autocorrelation_map(self, fn=None):
        ds = self.autocorrelation(self.field)
        if fn is not None:  ds.to_netcdf(fn)
        return ds
    
