import xarray as xr

from paths import path_samoc


class SST_index_analysis():
    """
    collection of analysis and plotting 
    """
    
    def __init__(self, index, run):
        # load data
        assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'had']
        assert index in ['AMO', 'SOM', 'TPI']
        
        self.run = run
        self.index = index
        self.data = xr.open_dataarray(f'{path_samoc}')
        self.raw_data = xr.open_dataarray(f'{path_samoc}')
        
    def spectrum(self):
       
    
    def fit_ar1(self):
        alpha = np.corrcoef(self.data.shift({time:1}),self.data)[0,1]
        
        
    def mc_confidence(self):
        
    
    def plot_regression_map(self):
        
    
    def plot_spectrum(self):
        
        
    def plot_timeseries_spectrum(self):
        """ plots both time series and spectrum side by side 
        spectrum plot includes null hypothesis confidence intervals
        """
        
        
class SST_index_comparison():
    
    def __init__(self, index):
        self.ctrl = xr.open_dataarray(f'{path_samoc}/SST/{index}_ctrl.nc')
        self.rcp  = xr.open_dataarray(f'{path_samoc}/SST/{index}_rcp.nc' )
        self.lpd  = xr.open_dataarray(f'{path_samoc}/SST/{index}_lpd.nc' )
        self.lpi  = xr.open_dataarray(f'{path_samoc}/SST/{index}_lpi.nc' )
        self.had  = xr.open_dataarray(f'{path_samoc}/SST/{index}_had.nc' )
        
    def plot_all_spectras(self):
        