import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import mtspec
from pandas.tools.plotting import autocorrelation_plot
from itertools import combinations

from maps import regr_map
from paths import path_samoc, path_results
from constants import A_earth
from analysis import TimeSeriesAnalysis, FieldAnalysis
# from xr_regression import lag_linregress_3D


class IndexAnalysis(TimeSeriesAnalysis):
    """
    collection of analysis and plotting functions
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
        """ loads final indices
        filtered and detrended with target region regression approach (Steinman et al. (2015))
        """
        def fn(run):
            return f'{path_samoc}/SST/{self.index}_{run}.nc'
        self.ctrl = xr.open_dataarray(fn('ctrl'), decode_times=False)
        self.rcp  = xr.open_dataarray(fn('rcp' ), decode_times=False)
        self.lpd  = xr.open_dataarray(fn('lpd' ), decode_times=False)
        self.lpi  = xr.open_dataarray(fn('lpi' ), decode_times=False)
        self.had  = xr.open_dataarray(fn('had' ), decode_times=False)
        self.all_indices = {'ctrl':self.ctrl, 'rcp':self.rcp, 'lpd':self.lpd, 'lpi':self.lpi, 'had':self.had}
        
    
    def load_raw_indices(self):
        """ loads raw, non-filtered indices
        detrended with target region regression approach (Steinman et al. (2015))
        """
        def fn_raw(run):
            return f'{path_samoc}/SST/{self.index}_raw_{run}.nc'
        self.ctrl_raw = xr.open_dataarray(fn_raw('ctrl'), decode_times=False)
        self.rcp_raw  = xr.open_dataarray(fn_raw('rcp' ), decode_times=False)
        self.lpd_raw  = xr.open_dataarray(fn_raw('lpd' ), decode_times=False)
        self.lpi_raw  = xr.open_dataarray(fn_raw('lpi' ), decode_times=False)
        if self.index in ['AMO', 'SOM']:
            self.had_raw  = xr.open_dataarray(f'{path_samoc}/SST/{self.index}_{self.index}_dt_raw_had.nc' , decode_times=False)
        elif self.index=='TPI':
            self.had_raw  = xr.open_dataarray(fn_raw('had') , decode_times=False)
        self.all_raw_indices = {'ctrl':self.ctrl_raw, 'rcp':self.rcp_raw, 'lpd':self.lpd_raw, 'lpi':self.lpi_raw, 'had':self.had_raw}
        
    
    def load_GMST_detr_indices(self):
        """ loads final indices
        filtered and detrended with scaled GMST approach (Kajtar et al. (2019))
        """
        def fn(run):
            return f'{path_samoc}/SST/{self.index}_GMST_dt_{run}.nc'
        self.ctrl_GMST = xr.open_dataarray(fn('ctrl'), decode_times=False)
        self.rcp_GMST  = xr.open_dataarray(fn('rcp' ), decode_times=False)
        self.lpd_GMST  = xr.open_dataarray(fn('lpd' ), decode_times=False)
        self.lpi_GMST  = xr.open_dataarray(fn('lpi' ), decode_times=False)
        self.had_GMST  = xr.open_dataarray(fn('had' ), decode_times=False)
        self.all_GMST_indices = {'ctrl':self.ctrl_GMST, 
                                 'rcp':self.rcp_GMST,
                                 'lpd':self.lpd_GMST,
                                 'lpi':self.lpi_GMST,
                                 'had':self.had_GMST}
        
        
    def load_yrly_SST(self):
        """ loads annual SST files """
        def fn(run):
            return f'{path_samoc}/SST/SST_yrly_{run}.nc'
        self.SST_ctrl = xr.open_dataarray(fn('ctrl'), decode_times=False)
        self.SST_rcp  = xr.open_dataarray(fn('rcp' ), decode_times=False)
        self.SST_lpd  = xr.open_dataarray(fn('lpd' ), decode_times=False)
        self.SST_lpi  = xr.open_dataarray(fn('lpi' ), decode_times=False)
        self.SST_had  = xr.open_dataarray(fn('had' ), decode_times=False)
        self.all_SSTs = {'ctrl':self.SST_ctrl,
                         'rcp':self.SST_rcp,
                         'lpd':self.SST_lpd,
                         'lpi':self.SST_lpi, 
                         'had':self.SST_had}
    
    
    def load_regression_files(self):
        """ loads regression files """
        def fn(run):
            return f'{path_samoc}/SST/{self.index}_regr_{run}.nc'
        self.regr_ctrl = xr.open_dataset(fn('ctrl'), decode_times=False)
        self.regr_rcp  = xr.open_dataset(fn('rcp' ), decode_times=False)
        self.regr_lpd  = xr.open_dataset(fn('lpd' ), decode_times=False)
        self.regr_lpi  = xr.open_dataset(fn('lpi' ), decode_times=False)
        self.regr_had  = xr.open_dataset(fn('had' ), decode_times=False)
        self.all_regrs = {'ctrl':self.regr_ctrl,
                          'rcp':self.regr_rcp,
                          'lpd':self.regr_lpd,
                          'lpi':self.regr_lpi,
                          'had':self.regr_had}
        
        
    def load_SST_autocorrelation_files(self):
        """ loads autocorrelation files """
        def fn(run):
            return f'{path_samoc}/SST/SST_autocorrelation_{run}.nc'
        self.autocorr_ctrl = xr.open_dataarray(fn('ctrl'), decode_times=False)
        self.autocorr_rcp  = xr.open_dataarray(fn('rcp' ), decode_times=False)
        self.autocorr_lpd  = xr.open_dataarray(fn('lpd' ), decode_times=False)
        self.autocorr_lpi  = xr.open_dataarray(fn('lpi' ), decode_times=False)
        self.autocorr_had  = xr.open_dataarray(fn('had' ), decode_times=False)
        self.all_autocorrs = {'ctrl':self.autocorr_ctrl, 
                              'rcp':self.autocorr_rcp,
                              'lpd':self.autocorr_lpd,
                              'lpi':self.autocorr_lpi,
                              'had':self.autocorr_had}

        
    def make_regression_files(self):
        """ generate regression files """
        self.load_SST_autocorrelation_files()
        for i, run in enumerate(self.all_indices.keys()):
            print(run)
            if run in ['ctrl', 'rcp', 'lpd', 'lpi']:  continue
            ds = self.lag_linregress(x=self.all_SSTs[run],
                                     y=self.all_indices[run],
                                     autocorrelation=self.all_autocorrs[run])
            ds.to_netcdf(f'{path_samoc}/SST/{self.index}_regr_{run}.nc')
        
        
    def plot_all_indices(self):
        """ plot time series of all indices """
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
        """ plot spectral estimation """
        fig, ax = plt.subplots(1, 1, figsize=(8,5))
        ax.tick_params(labelsize=14)
        ax.set_yscale('log')
        for i, run in enumerate(self.all_indices.keys()):
            spec, freq, jackknife, _, _ = mtspec.mtspec(
                data=self.all_indices[run].values, delta=1., time_bandwidth=4,
                number_of_tapers=5, statistics=True)
            ax.plot(freq, spec, color=f'C{i}', label=run.upper())
            ax.fill_between(freq, jackknife[:, 0], jackknife[:, 1], color=f'C{i}', alpha=0.1)
        ax.legend(fontsize=14, loc=1, frameon=False)
        ax.set_ylim((1E-6, 1E1))
        ax.set_xlim((0, 0.1))
        ax.set_xlabel('Period [yr]', fontsize=14)
        ax.set_xlabel(r'Frequency [yr$^{-1}$]', fontsize=14)
        ax.set_ylabel(f'{self.index} Power Spectral Density', fontsize=14)
        plt.tight_layout()
        plt.savefig(f'{path_results}/SST/{self.index}_all_spectra')
        
    
    def plot_spectrum_ar1(self, run):
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
        
        
    def plot_all_autocorrelations(self, n=60):
        fig, ax = plt.subplots(1, 1, figsize=(8,5))
        ax.tick_params(labelsize=14)
        ax.axhline(0, c='k', lw=.5)
        for i, run in enumerate(self.all_indices.keys()):
            tsa = TimeSeriesAnalysis(self.all_indices[run].values)
            acs = tsa.autocorrelation(n=n)
            plt.plot(np.arange(0,n+1), acs, label=f'{run.upper()}')
        plt.legend(fontsize=14, frameon=False)
        ax.set_xlim((0,n))
        ax.set_ylim((-1.1,1.1))
        ax.set_xlabel(r'lag [yr]', fontsize=14)
        ax.set_ylabel(f'{self.index} Autocorrelation', fontsize=14)
        plt.tight_layout()
        plt.savefig(f'{path_results}/SST/{self.index}_all_autocorrelations')
    
    
    def plot_regression_map(self, run):
        self.load_regression_files()
        regr_map(ds=self.all_regrs[run] , index=self.index, run=run , fn=None)
    
    
        
        
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
        self.load_TOA()
        self.load_SST_indices()
        
        
    def load_TOA(self):
        self.ctrl['FTNT'] = xr.open_dataarray(f'{path_samoc}/TOA/TOM_ctrl.nc', decode_times=False)/A_earth
        self.rcp ['FTNT'] = xr.open_dataarray(f'{path_samoc}/TOA/TOM_rcp.nc' , decode_times=False)/A_earth
        self.lpd ['FTNT'] = xr.open_dataarray(f'{path_samoc}/TOA/TOM_lpd.nc' , decode_times=False)/A_earth
        self.lpi ['FTNT'] = xr.open_dataarray(f'{path_samoc}/TOA/TOM_lpi.nc' , decode_times=False)/A_earth
        
        self.ctrl['FTNT_dt'] = self.ctrl['FTNT'] - self.lintrend( self.ctrl['FTNT'])
        self.rcp ['FTNT_dt'] = self.rcp ['FTNT'] - self.quadtrend(self.rcp ['FTNT'])
        self.lpd ['FTNT_dt'] = self.lpd ['FTNT'] - self.lintrend( self.lpd ['FTNT'])
        self.lpi ['FTNT_dt'] = self.lpi ['FTNT'] - self.lintrend( self.lpi ['FTNT'])
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
        timeseries = np.array([[self.ctrl['FTNT'].values, self.rcp['FTNT'].values,
                                self.lpd['FTNT'].values, self.lpi['FTNT'].values],
                               [self.ctrl['FTNT_dt'].values, self.rcp['FTNT_dt'].values,
                                self.lpd['FTNT_dt'].values, self.lpi['FTNT_dt'].values]])
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
    
    
class FieldSynthesis(FieldAnalysis):
    """ compares xr fields """
    def __init__(self, run):
        self.run = run
        
#     def spatial_