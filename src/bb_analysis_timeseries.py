"""
contains 3 classes:
"""

import numpy as np
import xesmf as xe
import mtspec
import xarray as xr

import scipy.stats as stats
import matplotlib.pyplot as plt
from statsmodels.tsa.arima_process import ArmaProcess
from statsmodels.stats.weightstats import DescrStatsW

from grid import generate_lats_lons
from regions import boolean_mask
from filters import lowpass, chebychev
from xr_DataArrays import dll_coords_names, xr_AREA

from ba_analysis_dataarrays import AnalyzeDataArray

    
    
class AnalyzeTimeSeries(AnalyzeDataArray):
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
#         assert type(self.ts)==np.ndarray
        assert 'time' in ts.coords
        self.len = len(self.ts)
        
        
    def rolling_trends(self, window=11):
        """ returns trends over a rolling window """
        assert type(window)==int
        da = xr.DataArray(data=np.full((self.len), np.nan),
                          coords={'time': self.ts.time}, dims=('time'))
        for t in range(self.len-window):
            da[int(np.floor(window/2))+t] = np.polyfit(np.arange(window), self.ts[t:t+window], 1)[0]
        return da
    
    
    def spectrum(self, data=None, filter_type=None, filter_cutoff=None):
        """ multitaper spectrum """
        if data is None:  data = np.array(self.ts)
        assert type(data)==np.ndarray
        if filter_type is not None:
            assert filter_type in ['lowpass', 'chebychev']
            assert type(filter_cutoff)==int
            assert filter_cutoff>1
            n = int(filter_cutoff/2)+1  # datapoints to remove from either end due to filter edge effects
            if filter_type=='lowpass':
                data = lowpass(data, filter_cutoff)[n:-n]
            elif filter_type=='chebychev':
                data = chebychev(data, filter_cutoff)[n:-n]
        
        spec, freq, jackknife, _, _ = mtspec.mtspec(
                data=data, delta=1., time_bandwidth=4,
                number_of_tapers=5, statistics=True)
        return (spec, freq, jackknife)
    
    
    def autocorrelation(self, data=None, n=1):
        """ calculates the first n lag autocorrelation coefficient """
        if data is None:
            assert n<self.len
            data = self.ts
        
        n += 1  # zeroth element is lag-0 autocorrelation = 1
        acs = np.ones((n))
        for i in np.arange(1,n):
            acs[i] = np.corrcoef(data[:-i]-data[:-i].mean(), data[i:]-data[i:].mean())[0,1]
        return acs
    
        
    def mc_ar1(self, n=1000):
        """ Monte-Carlo AR(1) processes """
        phi = self.autocorrelation(n=1)[1]
        AR_object = ArmaProcess(np.array([1, -phi]), np.array([1]))
#         sigma_eps = np.sqrt(np.var(self.ts)*(1-phi**2))
#         AR_object = ArmaProcess(np.array([1, -phi]), np.array([1, sigma_eps]))
        mc = np.zeros((n, self.len))
        for i in range(n):
            mc[i,:] = AR_object.generate_sample(nsample=self.len)
            mc[i,:] *= np.std(self.ts.values)/np.std(mc[i,:])
        return mc
    
        
    def mc_ar1_spectrum(self, N=1000, filter_type=None, filter_cutoff=None):
        """ calculates the MC avg spectrum and the 95% confidence interval """
        mc = self.mc_ar1()
        if filter_type is not None:
            assert filter_type in ['lowpass', 'chebychev']
            assert type(filter_cutoff)==int
            assert filter_cutoff>1
            n = int(filter_cutoff/2)+1  # datapoints to remove from either end due to filter edge effects
            mc = mc[:,n:-n]
            if filter_type=='lowpass':      mc = lowpass(mc.T, filter_cutoff).T
            elif filter_type=='chebychev':  mc = chebychev(mc.T, filter_cutoff).T

        mc_spectra = np.zeros((N, int(len(mc[0,:])/2)+1))#int(self.len/2)+1))
        for i in range(N):
            (mc_spectra[i,:], freq, jk) = self.spectrum(data=mc[i,:])
        mc_spectrum = np.zeros((4, int(len(mc[0,:])/2)+1))#int(self.len/2)+1))
        mc_spectrum[0,:] = np.median(mc_spectra, axis=0)
        mc_spectrum[1,:] = freq
        mc_spectrum[2,:] = np.percentile(mc_spectra,  5, axis=0)
        mc_spectrum[3,:] = np.percentile(mc_spectra, 95, axis=0)
        return mc_spectrum
    
    
    @staticmethod
    def test_homoscedasticity(X, Y):
        X1, Y1 = X, Y
        if len(X)>150:  X1 = X[-150:]
        if len(Y)>150:  Y1 = Y[-150:]
        print(f'{stats.levene(X, Y)[1]:4.2e}, {stats.levene(X1, Y1)[1]:4.2e}')
        
        
    def plot_spectrum_ar1(self, data=None):
        """ plots spectrum of single time series + AR(1) spectrum
        AR(1) spectrum includes uncertainties from MC simulation
        for the filtered SST indices, the AR(1) process is first fitted to the annual data
        """
        self.load_raw_indices()
        tsa = TimeSeriesAnalysis(self.all_raw_indices[run].values)
        ft, fc = 'lowpass', 13
        spectrum = tsa.spectrum(filter_type=ft, filter_cutoff=fc)
        mc_spectrum = tsa.mc_ar1_spectrum(filter_type=ft, filter_cutoff=fc)
        
        fig, ax = plt.subplots(1, 1, figsize=(8,5))
        ax.tick_params(labelsize=14)
        ax.set_yscale('log')        
        L2 = ax.fill_between(mc_spectrum[1,:], mc_spectrum[2,:],  mc_spectrum[3,:],
                        color='C1', alpha=.3, label='5-95% C.I.')
        L1, = ax.plot(mc_spectrum[1,:], mc_spectrum[0,:], c='C1', label=f'MC AR(1)')     
        L4 = ax.fill_between(spectrum[1], spectrum[2][:, 0], spectrum[2][:, 1],
                       color='C0', alpha=0.3, label=f'{run.upper()} jackknife estimator')
        L3, = ax.plot(spectrum[1], spectrum[0], c='C0', label=f'{run.upper()} spectrum')
        leg1 = plt.legend(handles=[L1, L2], fontsize=14, frameon=False, loc=3)
        ax.legend(handles=[L3,L4], fontsize=14, frameon=False, loc=1)
        ax.add_artist(leg1)
        ymax = 1e1
        if any([mc_spectrum[3,:].max()>ymax, spectrum[2][:, 1].max()>ymax]):
            ymax = max([mc_spectrum[3,:].max(), spectrum[2][:, 1].max()])
        ax.set_ylim((1E-6, ymax))
        ax.set_xlim((0, 0.1))
        ax.set_xlabel(r'Frequency [yr$^{-1}$]', fontsize=14)
        ax.set_ylabel(f'{self.index} Power Spectral Density', fontsize=14)
        plt.tight_layout()
        plt.savefig(f'{path_results}/SST/{self.index}_AR1_spectrum_{run}')
        