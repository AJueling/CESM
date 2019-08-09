import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import mtspec
from pandas.tools.plotting import autocorrelation_plot
from itertools import combinations

from maps import regr_map
from paths import path_samoc, path_results
from constants import A_earth
from bb_analysis_timeseries import TimeSeriesAnalysis
from bc_analysis_fields import FieldAnalysis

# =============================================================================
# =============================================================================

        
class TimeSeriesSynthesis(TimeSeriesAnalysis):
    """ combines different time series """
    def __init__(self):
        self.ctrl = {'name':'ctrl'}
        self.rcp  = {'name':'rcp'}
        self.lpd  = {'name':'lpd'}
        self.lpi  = {'name':'lpi'}
        self.had  = {'name':'had'}
        
        self.runs = [self.ctrl, self.rcp, self.lpd, self.lpi, self.had]
        
        self.load_GMST()
#         self.load_TOA()
        self.load_SST_indices()
        
        
    def load_TOA(self):
        def fn(run):
            return f'{path_samoc}/TOA/TOM_{run}.nc'
        self.ctrl['FTNT'] = xr.open_dataarray(fn('ctrl'), decode_times=False)\
                            /A_earth
        self.rcp ['FTNT'] = xr.open_dataarray(fn('rcp') , decode_times=False)\
                            /A_earth
        self.lpd ['FTNT'] = xr.open_dataarray(fn('lpd') , decode_times=False)\
                            /A_earth
        self.lpi ['FTNT'] = xr.open_dataarray(fn('lpi') , decode_times=False)\
                            /A_earth
        
        self.ctrl['FTNT_dt'] = self.ctrl['FTNT']\
                               - self.lintrend( self.ctrl['FTNT'])
        self.rcp ['FTNT_dt'] = self.rcp ['FTNT']\
                               - self.quadtrend(self.rcp ['FTNT'])
        self.lpd ['FTNT_dt'] = self.lpd ['FTNT']\
                               - self.lintrend( self.lpd ['FTNT'])
        self.lpi ['FTNT_dt'] = self.lpi ['FTNT']\
                               - self.lintrend( self.lpi ['FTNT'])
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
        
        
#     @staticmethod
    def print_pairwise_homoscedasticity(self, fields):
        assert type(fields)==tuple
        print('p-value\nfull time series, last 150 years')
        for c in combinations(fields, 2):
            self.test_homoscedasticity(c[0], c[1])
        return
            
    def print_all_autocorrelations(self, timeseries):
        for ts in timeseries:
            print(self.autocorrelation(ts))
        return
    
    def model_GMST_modes(self):
        """ linear models of GMST as a function of SST indices """
        return
    
    def model_GMST_TOA_modes(self):
        """ linear models of GMST as function of SST modes and TOA """
        return
    
    def model_GMST_TOA_OHU(self):
        """ GMST as a function of TOA and OHU """
        return
    
    
    @staticmethod
    def plot_all_timeseries(times, timeseries, ylabels, xlim, fn=None):
        (m,n) = np.shape(timeseries)
        
        f, ax, = plt.subplots(m,1,figsize=(12,m*3), sharex=True)
        for i in range(m):
            ax[i].tick_params(labelsize=14)
            ax[i].axhline(0, c='k', lw=.5)
        # ax.axhline(0,c='k', lw=.5)

        time_had = np.arange(2350,2519)

        for i in range(m):
            for j in range(n):
                time = times[i,j]
                data = timeseries[i,j]
                ax[i].plot(time, data, c=f'C{j}')
                ax[i].set_ylabel(ylabels[i], fontsize=16)
                
        ax[-1].set_xlabel('time [years]', fontsize=16)

        ax[-1].set_xticks(np.arange(0,2800,100))
        ax[-1].set_xlim(xlim)
        f.align_ylabels()
        f.tight_layout()
        if fn is not None:  plt.savefig(fn)
            
            
    def plot_all_spectra(self, timeseries, fn=None):
        f, ax = plt.subplots(1,1, figsize=(8,5))
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.tick_params(labelsize=14)
        
        for ts in timeseries:
            spectrum = self.spectrum(data=ts.values)
            ax.plot(spectrum[1], spectrum[0])
        ax.set_xlabel(r'frequency [yr$^{-1}$]', fontsize=14)
        ax.set_ylabel('power spectral density', fontsize=14)
        plt.tight_layout()
        
        
    def plot_all_TOA(self):
        self.load_TOA()
        time_ctrl = np.arange(1950,2150)
        time_rcp  = np.arange(2200,2300)
        time_lpd  = np.arange(1503,1851)
        time_lpi  = np.arange(1276,1482)
        times = np.array([time_ctrl, time_rcp, time_lpd, time_lpi])
        times = np.array([times, times])
        timeseries = np.array([[self.ctrl['FTNT'].values,
                                self.rcp['FTNT'].values,
                                self.lpd['FTNT'].values,
                                self.lpi['FTNT'].values],
                               [self.ctrl['FTNT_dt'].values,
                                self.rcp['FTNT_dt'].values,
                                self.lpd['FTNT_dt'].values, 
                                self.lpi['FTNT_dt'].values]])
        ylabels = [r'TOA [W m$^{-2}$]', r'detrended TOA [W m$^{-2}$]']
        xlim = (1230,2350)
        fn = f'{path_results}/TOA/TOA_timeseries'
        self.plot_all_timeseries(times, timeseries, ylabels, xlim, fn=fn)
        
    def print_TOA_statistics(self):
        self.load_TOA()
        for i, ts in enumerate([self.ctrl['FTNT_dt'], self.rcp['FTNT_dt'],\
                                self.lpd['FTNT_dt'], self.lpi['FTNT_dt']]):
            print(i, TimeSeriesAnalysis(ts).autocorrelation(), ts.std().values)
        
        
    def plot_lead_lag(self):
        return
    
    
    def plot_variance(self):
        return