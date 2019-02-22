import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import mtspec
from pandas.tools.plotting import autocorrelation_plot

from maps import regr_map
from paths import path_samoc, path_results
from xr_regression import lag_linregress_3D

class time_series_analysis(object):
    
    def __init__(self, ts=None):
        """ """
        if ts=='GMST':
            self.ts = xr.open_dataarray('')
        elif ts in ['AMO', 'SOM', 'TPI']:
            self.ts = xr.open_dataarray('')
    
    def spectrum(self):
        return
    
    def fit_ar1(self):
        alpha = np.corrcoef(self.data.shift({time:1}),self.data)[0,1]
        return
        
    def mc_confidence(self):
        return
        
        
    
    def plot_autocorrelation(self):
        return
    
    
    def plot_spectrum(self):
        return
        
    
    def plot_timeseries_spectrum(self):
        """ plots both time series and spectrum side by side 
        spectrum plot includes null hypothesis confidence intervals
        """
        return
    
    
    def plot_regression_map(self):
        return
    
    
    

class SST_index_analysis(time_series_analysis):
    """
    collection of analysis and plotting 
    """
    
    def __init__(self, index, run='all'):
        assert run in ['all', 'ctrl', 'rcp', 'lpd', 'lpi', 'had']
        assert index in ['AMO', 'SOM', 'TPI']
        
        self.index = index
        self.load_indices()
        self.load_GMST_detr_indices()
        self.load_yrly_SST()
        
        self.plot_time_offsets = [1850, 200, 1350, -1600, 2350]
        self.plot_texts = ['CTRL' ,'RCP' ,'pres. day low' ,'pre-ind. low' ,'HadISST']
        self.plot_text_loc = [1950, 2200,1500 ,0 ,2320]
        
    def create_SST_avg_index(self, run, detrend_signal):
        # copy SST_index function
        return
        
        
    def load_indices(self):
        """ loads final indices filtered and detrended with target region regression approach (Steinman et al. (2015))"""
        self.ctrl = xr.open_dataarray(f'{path_samoc}/SST/{self.index}_ctrl.nc', decode_times=False)
        self.rcp  = xr.open_dataarray(f'{path_samoc}/SST/{self.index}_rcp.nc' , decode_times=False)
        self.lpd  = xr.open_dataarray(f'{path_samoc}/SST/{self.index}_lpd.nc' , decode_times=False)
        self.lpi  = xr.open_dataarray(f'{path_samoc}/SST/{self.index}_lpi.nc' , decode_times=False)
        self.had  = xr.open_dataarray(f'{path_samoc}/SST/{self.index}_had.nc' , decode_times=False)
        self.all_indices = {'ctrl':self.ctrl, 'rcp':self.rcp, 'lpd':self.lpd, 'lpi':self.lpi, 'had':self.had}
        
    
    def load_GMST_detr_indices(self):
        """ loads final indices filtered and detrended with scaled GMST approach (Kajtar et al. (2019))"""
        self.ctrl_GMST = xr.open_dataarray(f'{path_samoc}/SST/{self.index}_GMST_dt_ctrl.nc', decode_times=False)
        self.rcp_GMST  = xr.open_dataarray(f'{path_samoc}/SST/{self.index}_GMST_dt_rcp.nc' , decode_times=False)
        self.lpd_GMST  = xr.open_dataarray(f'{path_samoc}/SST/{self.index}_GMST_dt_lpd.nc' , decode_times=False)
        self.lpi_GMST  = xr.open_dataarray(f'{path_samoc}/SST/{self.index}_GMST_dt_lpi.nc' , decode_times=False)
        self.had_GMST  = xr.open_dataarray(f'{path_samoc}/SST/{self.index}_GMST_dt_had.nc' , decode_times=False)
        self.all_GMST_indices = {'ctrl':self.ctrl_GMST, 'rcp':self.rcp_GMST, 'lpd':self.lpd_GMST, 'lpi':self.lpi_GMST, 'had':self.had_GMST}
        
        
    def load_yrly_SST(self):
        self.SST_ctrl = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_ctrl.nc', decode_times=False)
        self.SST_rcp  = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_rcp.nc' , decode_times=False)
        self.SST_lpd  = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_lpd.nc' , decode_times=False)
        self.SST_lpi  = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_lpi.nc' , decode_times=False)
        self.SST_had  = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_had.nc' , decode_times=False)
        self.all_SSTs = {'ctrl':self.SST_ctrl, 'rcp':self.SST_rcp, 'lpd':self.SST_lpd, 'lpi':self.SST_lpi, 'had':self.SST_had}
    
    
    def load_regression_files(self):
        self.regr_ctrl = xr.open_dataset(f'{path_samoc}/SST/{self.index}_regr_ctrl.nc', decode_times=False)
        self.regr_rcp  = xr.open_dataset(f'{path_samoc}/SST/{self.index}_regr_rcp.nc' , decode_times=False)
        self.regr_lpd  = xr.open_dataset(f'{path_samoc}/SST/{self.index}_regr_lpd.nc' , decode_times=False)
        self.regr_lpi  = xr.open_dataset(f'{path_samoc}/SST/{self.index}_regr_lpi.nc' , decode_times=False)
        self.regr_had  = xr.open_dataset(f'{path_samoc}/SST/{self.index}_regr_had.nc' , decode_times=False)
        self.all_regrs = {'ctrl':self.regr_ctrl, 'rcp':self.regr_rcp, 'lpd':self.regr_lpd, 'lpi':self.regr_lpi, 'had':self.regr_had}

        
    def make_regression_files(self):
        for i, run in enumerate(self.all_indices.keys()):
            print(run)
            ds = lag_linregress_3D(x=self.all_SSTs[run], y=self.all_indices[run])
            ds.to_netcdf(f'{path_samoc}/SST/{self.index}_regr_{run}.nc')
        
        
    def plot_all_indices(self):
        assert self.all_indices and self.all_GMST_indices
        f, ax = plt.subplots(1,1,figsize=(12,5))
        ax.tick_params(labelsize=14)
        ax.axhline(0, c='k', lw=.5)
        maxv = 0
        for i, run in enumerate(self.all_indices.keys()):
            if max(self.all_indices[run])>maxv:
                maxv = max(self.all_indices[run])
            dt = self.plot_time_offsets[i]
            L1, = ax.plot(self.all_indices[run].time/365+dt, self.all_indices[run],
                          c=f'C{i}', ls='-', lw=1, label='SST detrended with target region regression')
            L2, = ax.plot(self.all_GMST_indices[run].time/365+dt, self.all_GMST_indices[run],
                          c=f'C{i}', ls='--', lw=1, label='SST detrended with scaled GMST') 
        for i in range(5):
            text = self.plot_texts[i]
            tloc = self.plot_text_loc[i]
            ax.text(tloc, maxv, text, fontsize=16, color=f'C{i}')
        ax.legend(handles=[L1,L2], loc=8, ncol=3, fontsize=14, frameon=False)
        ax.set_ylabel(f'{self.index} index [K]', fontsize=16)
        ax.set_xlabel('time [years]', fontsize=16)
        ax.set_xticks(np.arange(0,2700,200))
        ax.set_xlim((-50,2550))
        f.tight_layout()
        f.savefig(f'{path_results}/SST/{self.index}_index_overview.png')    
    
    
    def plot_all_spectra(self):
        """plot spectral estimation """
        fig, ax = plt.subplots(1, 1, figsize=(8,5))
        ax.tick_params(labelsize=14)
        ax.set_yscale('log')
#         ax.set_xscale('log')
        for i, run in enumerate(self.all_indices.keys()):
            spec, freq, jackknife, _, _ = mtspec.mtspec(
                data=self.all_indices[run].values, delta=1., time_bandwidth=4,
                number_of_tapers=5, statistics=True)
            ax.plot(freq, spec, color=f'C{i}', label=run.upper())
            ax.fill_between(freq, jackknife[:, 0], jackknife[:, 1], color=f'C{i}', alpha=0.1)
        ax.legend(fontsize=14, loc=1, frameon=False)
        ax.set_ylim((1E-6, 1E1))
#         ax.set_xlim((5, 2e3))
        ax.set_xlim((0, 0.1))
        ax.set_xlabel('Period [yr]', fontsize=14)
        ax.set_xlabel(r'Frequency [yr$^{-1}$]', fontsize=14)
        ax.set_ylabel(f'{self.index} Power Spectral Density', fontsize=14)
        plt.tight_layout()
        
        
    def plot_all_autocorrelations(self):
#         fig, ax = plt.subplots(1, 1, figsize=(8,5))
#         ax.tick_params(labelsize=14)
        for i, run in enumerate(self.all_indices.keys()):
#             fig, ax = plt.subplots(1, 1, figsize=(8,5))
#             ax.tick_params(labelsize=14)
            autocorrelation_plot(self.all_indices[run].values, figsize=(8,5))
#             ax.set_xlim((0,100))
        
        return
    
    
    def plot_all_regression_maps(self):
        self.load_regression_files()
#         for i, run in enumerate(self.all_indices.keys()):
        regr_map(ds=self.regr_had, index=self.index, run='had', fn=None)
    
    
        