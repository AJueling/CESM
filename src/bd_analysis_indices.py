import glob
import os
import xarray as xr
import matplotlib.pyplot as plt

from paths import path_samoc
from regions import boolean_mask, TPI_masks, mask_box_in_region,\
                    bll_AMO, bll_SOM, bll_TPI1, bll_TPI2, bll_TPI3
from filters import chebychev, lowpass
from timeseries import IterateOutputCESM
from xr_regression import xr_lintrend
from xr_DataArrays import xr_AREA

from ab_derivation_SST import DeriveSST
from ba_analysis_dataarrays import AnalyzeDataArray
from bc_analysis_fields import AnalyzeField

class AnalyzeIndex(object):
    """ calculating SST indices """
    def __init__(self):
        return
    
    
    def SST_area_average(self, xa_SST, AREA, AREA_index, MASK, dims=('nlat', 'nlon'), index_loc=None):
        """ calculates the average SST over an area, possibly as a time series """
        assert type(xa_SST)==xr.core.dataarray.DataArray
        assert type(AREA)==xr.core.dataarray.DataArray
        print(f'calculating average of SST in area {index_loc}')
        if type(index_loc)==dict:
            index = (xa_SST*AREA).where(MASK).sel(index_loc).sum(dim=dims)/AREA_index
        elif index_loc==None:
            index = (xa_SST*AREA).where(MASK).sum(dim=dims)/AREA_index
        else:
            print('kwarg `index_loc` is not given properly.')
        return index
    
    
    def SST_index(self, index, run, detrend_signal='GMST', time_slice='full'):
        """ calcalates SST time series from yearly detrended SST dataset """
        assert index in ['AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']
        assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'had']
        assert detrend_signal in ['GMST', 'AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']

        print(index, run)

        if run in ['ctrl', 'rcp']:
            domain = 'ocn'  #check this
            dims = ('nlat', 'nlon')
        elif run in ['lpd', 'lpi']:
            domain = 'ocn_low'
            dims = ('nlat', 'nlon')
        elif run=='had':
            domain = 'ocn_had'
            dims = ('latitude', 'longitude')

        blats, blons, mask_nr = self.bounding_lats_lons(index)
        MASK = mask_box_in_region(domain=domain, mask_nr=mask_nr,\
                                  bounding_lats=blats, bounding_lons=blons)
        AREA = xr_AREA(domain=domain).where(MASK)
        index_area = AREA.sum()

        if detrend_signal=='GMST' or run=='had' and detrend_signal in ['AMO', 'SOM']:
            print(f'underlying SST field: detrended with {detrend_signal}, no filtering')
            if detrend_signal=='GMST':
                print('GMST(t) signal scaled at each grid point')
                dt = self.detrend_string(run)
            else:
                print(f'{detrend_signal}(t) removed from all SST gridpoints without scaling\n')
            detr = f'_{detrend_signal}_dt'
                
            fn = f'{path_samoc}/SST/SST_{detrend_signal}_{dt}_yrly_{run}.nc'
            assert os.path.exists(fn)
            SST_yrly = xr.open_dataarray(fn).where(MASK)
            if time_slice!='full':
                SST_yrly = DeriveSST.select_time_slice(SST_yrly, time_slice)
            print('after loading of data')
            
        else:  # run=='had' and detrend_signal!='GMST'
            print('underlying SST field: no detrending, no filtering')
            if detrend_signal in ['AMO', 'SOM']:
                print(f'{detrend_signal} must subsequently be detrended with polynomial\n')
            else:
                print(f'{detrend_signal} must not be detrended since forcing signal compensated in TPI\n')
            SST_yrly = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_{run}.nc').where(MASK)
            detr = ''

        SSTindex = self.SST_area_average(xa_SST=SST_yrly, AREA=AREA,\
                                         AREA_index=index_area, MASK=MASK, dims=dims)
        ts = self.time_slice_string(time_slice)
        fn = f'{path_samoc}/SST/{index}_{detrend_signal}_{dt}_raw_{run}{ts}.nc'
        SSTindex.to_netcdf(fn)
        return SSTindex
    
    
    def derive_all_SST_avg_indices(self, run, tslice):
        """ generates all SST avg indices detrended with the GMST signal  for full time series """
        AMO  = self.SST_index('AMO' , run, detrend_signal='GMST', time_slice=tslice)
        SOM  = self.SST_index('SOM' , run, detrend_signal='GMST', time_slice=tslice)
        TPI1 = self.SST_index('TPI1', run, detrend_signal='GMST', time_slice=tslice)
        TPI2 = self.SST_index('TPI2', run, detrend_signal='GMST', time_slice=tslice)
        TPI3 = self.SST_index('TPI3', run, detrend_signal='GMST', time_slice=tslice)
        return
  
        
    def derive_final_SST_indices(self, run, tslice):
        """ processes raw indices: filtering and TPI summation """
        assert run in ['ctrl', 'lpd', 'had']
        assert tslice=='full' or type(tslice)==tuple
        
        ts = self.time_slice_string(tslice)
        dt = self.detrend_string(run)
            
        # AMO & SOM
        for i, idx in enumerate(['AMO', 'SOM']):
            fn = f'{path_samoc}/SST/{idx}_GMST_{dt}_raw_{run}.nc'
            fn_new = f'{path_samoc}/SST/{idx}_{run}{ts}.nc'
            da = xr.open_dataarray(fn, decode_times=False)
            da = DeriveSST.select_time_slice(da, tslice)
            lowpass(da, 13).to_netcdf(fn_new)

        # TPI
        fn1 = f'{path_samoc}/SST/TPI1_GMST_{dt}_raw_{run}.nc'
        TPI1 = xr.open_dataarray(fn1, decode_times=False)
        fn2 = f'{path_samoc}/SST/TPI2_GMST_{dt}_raw_{run}.nc'
        TPI2 = xr.open_dataarray(fn2, decode_times=False)
        fn3 = f'{path_samoc}/SST/TPI3_GMST_{dt}_raw_{run}.nc'
        TPI3 = xr.open_dataarray(fn3, decode_times=False)
        TPI = TPI2 - (TPI1+TPI3)/2
        TPI = DeriveSST.select_time_slice(TPI, tslice)
        lowpass(TPI, 13).to_netcdf(f'{path_samoc}/SST/TPI_{run}{ts}.nc')
        return
    
    
    def derive_yrly_autocorrelations(self, run, tslice):
        """ autocorrelation maps for detrended SST fields for significance tests """
        ts = self.time_slice_string(tslice)
        dt = self.detrend_string(run)
        
        fn = f'{path_samoc}/SST/SST_GMST_{dt}_yrly_{run}{ts}.nc'
        fn_new = f'{path_samoc}/SST/SST_autocorrelation_{run}{ts}.nc'

        da = xr.open_dataarray(fn, decode_times=False)
        FA = AnalyzeField(da)
        FA.make_autocorrelation_map(fn_new)
        return
    
    
    def make_yrly_regression_files(self, run, tslice):
        """ generate regression files """
        assert run in ['ctrl', 'lpd', 'had']
        
        ts = self.time_slice_string(tslice)
        dt = self.detrend_string(run)

        for idx in ['AMO', 'SOM', 'TPI']:
            fn_idx = f'{path_samoc}/SST/{idx}_{run}{ts}.nc'
            fn_SST = f'{path_samoc}/SST/SST_GMST_{dt}_yrly_{run}{ts}.nc'
            fn_acr = f'{path_samoc}/SST/SST_autocorrelation_{run}{ts}.nc'

            index = xr.open_dataarray(fn_idx, decode_times=False)
            SST_dt = xr.open_dataarray(fn_SST, decode_times=False)
            autocorr = xr.open_dataarray(fn_acr, decode_times=False)

            xA = AnalyzeDataArray()
            ds = xA.lag_linregress(x=index[7:-7],  # removing filter edge effects
                                   y=SST_dt[7:-7], 
                                   autocorrelation=autocorr,
                                   standardize=True,
                                  )
            ds.to_netcdf(f'{path_samoc}/SST/{idx}_regr_{run}{ts}.nc')
        print('success')
        return
    
    
    ## AUXILIARY FUNCTIONS
        
        
    def bounding_lats_lons(self, index):
        """ bounding latitudes and longitudes """
        if index=='AMO':
            (blats, blons) = bll_AMO
            mask_nr = 6
        elif index=='SOM':
            (blats, blons) = bll_SOM
            mask_nr = 0
        elif index=='TPI1':
            (blats, blons) = bll_TPI1
            mask_nr = 2
        elif index=='TPI2':
            (blats, blons) = bll_TPI3
            mask_nr = 2
        elif index=='TPI3':
            (blats, blons) = bll_TPI3
            mask_nr = 2
        return blats, blons, mask_nr
            
    
    def time_slice_string(self, tslice):
        """ string for time subset """
        if tslice=='full':         ts = ''
        elif type(tslice)==tuple:  ts = f'_{tslice[0]}_{tslice[1]}'
        else:                      raise ValueError()
        return ts
            
        
    def detrend_string(self, run):
        """ scaled linear/quadratic or two factor detrending string in filenames """
        if run=='had':
            dt = 'tfdt'  # two-factor detrending, or 'sfdt' single-factor detrending
        elif run in ['ctrl', 'lpd', 'rcp']:
            dt = 'sqdt'  # scaled quadratic detrending
        else:
            raise ValueError()
        return dt  
    
    
    
# =============================================================================
# GLOBAL MEAN TIME SERIES
# =============================================================================

# %%time
# # 29 min
# TAREA = xr_AREA('ocn')
# MASK_ocn = boolean_mask('ocn', 0)
# global_area = TAREA.where(MASK_ocn).sum()
# for run in ['ctrl', 'rcp']:
#     for i, (y,m,s) in enumerate(IterateOutputCESM(domain='ocn', run=run, tavg='monthly')):
#         SST = xr.open_dataset(s, decode_times=False).TEMP[0,0,:,:]
#         SST_gm = SST_index(xa_SST=SST, AREA=TAREA, index_loc=global_ocean, AREA_index=global_area, MASK=MASK_ocn, dims=('nlat', 'nlon'))
#         if m==1: print(y, SST_gm.item()) 
#         if i==0:  SST_new = SST_gm
#         else:     SST_new = xr.concat([SST_new, SST_gm], dim='time')
#     SST_new.to_netcdf(f'{path_samoc}/SST/SST_global_mean_monthly_{run}.nc')

# %%time
# # CTRL and RCP global means: ocn_rect
# # 7 min for both
# AREA_rect = xr_AREA('ocn_rect')
# MASK_rect = boolean_mask('ocn_rect', 0)
# global_area2 = AREA_rect.where(MASK_rect).sum()
# for run in ['ctrl', 'rcp']:
#     for i, (y,m,s) in enumerate(IterateOutputCESM(domain='ocn_rect', run=run, tavg='monthly')):
#         SST = xr.open_dataset(s, decode_times=False).TEMP[0,:,:]
#         SST_gm = SST_index(xa_SST=SST, AREA=AREA_rect, index_loc=gl_ocean_rect,
#                            AREA_index=global_area2, MASK=MASK_rect,
#                            dims=('t_lat', 't_lon'))
#         if m==1: print(y, SST_gm.item()) 
#         if i==0:  SST_new = SST_gm
#         else:     SST_new = xr.concat([SST_new, SST_gm], dim='time')
#     SST_new.to_netcdf(f'{path_samoc}/SST/SST_global_mean_monthly_rect_{run}.nc')

# %%time
# # CTRL and RCP 60S to 60N means: ocn_rect
# # 7:26 min for both
# AREA_rect = xr_AREA('ocn_rect').sel({'t_lat':slice(-60,60)})
# MASK_rect = boolean_mask('ocn_rect', 0).sel({'t_lat':slice(-60,60)})
# global_area = AREA_rect.where(MASK_rect).sum()
# for run in ['ctrl', 'rcp']:
#     for i, (y,m,s) in enumerate(IterateOutputCESM(domain='ocn_rect', run=run, tavg='monthly')):
#         SST = xr.open_dataset(s, decode_times=False).TEMP[0,:,:].sel({'t_lat':slice(-60,60)})
#         SST_gm = SST_index(xa_SST=SST, AREA=AREA_rect, index_loc=gl_ocean_rect,
#                            AREA_index=global_area, MASK=MASK_rect,
#                            dims=('t_lat', 't_lon'))
#         if m==1: print(y, SST_gm.item()) 
#         if i==0:  SST_new = SST_gm
#         else:     SST_new = xr.concat([SST_new, SST_gm], dim='time')
#     SST_new.to_netcdf(f'{path_samoc}/SST/SST_60S_60N_mean_monthly_rect_{run}.nc')




# =============================================================================
# 60S-60N MEAN TIME SERIES
# =============================================================================

# %%time
# # LPI and LPD global means: ocn_low
# # 1 hr for LPD, 19 min for LPI
# AREA_low = xr_AREA('ocn_low')
# MASK_low = boolean_mask('ocn_low', 0)
# global_area_low = AREA_low.where(MASK_low).sum()
# for run in ['lpi']:#, 'lpd']:
#     for i, (y,m,s) in enumerate(IterateOutputCESM(domain='ocn_low', run=run, tavg='monthly')):
#         SST = xr.open_dataset(s, decode_times=False).TEMP[0,0,:,:]
#         SST_gm = SST_index(xa_SST=SST, AREA=AREA_low, index_loc=gl_ocean_low, 
#                            AREA_index=global_area_low, MASK=MASK_low, 
#                            dims=('nlat', 'nlon'))
#         if m==1: print(y, SST_gm.item()) 
#         if i==0:  SST_new = SST_gm
#         else:     SST_new = xr.concat([SST_new, SST_gm], dim='time')
#     SST_new.to_netcdf(f'{path_samoc}/SST/SST_global_mean_monthly_{run}.nc')



# HadISST 60S-60N time series is generated in SST_obs.ipynb

# =============================================================================
# ISOLATE PACIFIC SST DATA
# =============================================================================








# ds = xr.open_dataset(file_HadISST, decode_times=False)
# run = 'had'
# for i, r in enumerate(['Pac_38S', 'Pac_Eq', 'Pac_20N']):
#     latS = [-38, 0, 20][i]
#     lonE = [300, 285, 255][i]
#     Pac_MASK = mask_box_in_region(domain='ocn_had', mask_nr=2, bounding_lats=(latS,68), bounding_lons=(110,lonE))
#     da = ds.sst.where(Pac_MASK, drop=False)
#     da.to_netcdf(f'{path_samoc}/SST/SST_monthly_{r}_{run}.nc')






# =============================================================================
# remove yearly files
# =============================================================================

# def cleanup_yearly_files():
# return