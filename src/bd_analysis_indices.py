import glob
import os
import xarray as xr
import matplotlib.pyplot as plt

from SST import SST_index
from paths import path_samoc
from timeseries import IterateOutputCESM
from xr_regression import xr_lintrend

yrly_ctrl_100_248 = (36531,90551)
yrly_ctrl_151_299 = (55146,109166)
yrly_lpd_268_416  = (97851,151871)
yrly_lpd_417_565  = (152236,206256)

from filters import chebychev, lowpass
from analysis import FieldAnalysis, xrAnalysis
from 2_1_analysis_dataarray import AnalyzeDataArray
from 2_3_analysis_field import AnalyzeField

class AnalyzeIndex(object):
    def __init__():
        return
    
    @staticmethod
    def derive_all_SST_avg_indices(run):
        """ generates all SST avg indices detrended with the GMST signal  for full time series """
        AMO  = SST_index('AMO' , run, detrend_signal='GMST')
        SOM  = SST_index('SOM' , run, detrend_signal='GMST')
        TPI1 = SST_index('TPI1', run, detrend_signal='GMST')
        TPI2 = SST_index('TPI2', run, detrend_signal='GMST')
        TPI3 = SST_index('TPI3', run, detrend_signal='GMST')
    
        if run=='ctrl':  slices = [yrly_ctrl_100_248, yrly_ctrl_151_299]
        elif run=='lpd': slices = [yrly_lpd_268_416, yrly_lpd_417_565]
            
        if run in ['ctrl', 'lpd']:
            for tslice in slices:
                AMO  = SST_index('AMO' , run, detrend_signal='GMST', time_slice=tslice)
                SOM  = SST_index('SOM' , run, detrend_signal='GMST', time_slice=tslice)
                TPI1 = SST_index('TPI1', run, detrend_signal='GMST', time_slice=tslice)
                TPI2 = SST_index('TPI2', run, detrend_signal='GMST', time_slice=tslice)
                TPI3 = SST_index('TPI3', run, detrend_signal='GMST', time_slice=tslice)
                
    @staticmethod
    def derive_final_SST_indices(run):
        """"""
        plt.figure()
        tslices = ['']
        if run=='ctrl':  tslices.extend(['_100_248', '_151_299'])
        elif run=='lpd':  tslices.extend(['_268_416', '_417_565'])
        
        for j, tslice in enumerate(tslices):
            ls= ['-', '--', '-.'][j]
            for i, idx in enumerate(['AMO', 'SOM']):    
                fn = f'{path_samoc}/SST/{idx}_GMST_dt_raw{tslice}_{run}.nc'
                fn_new = f'{path_samoc}/SST/{idx}{tslice}_{run}.nc'
                da = xr.open_dataarray(fn, decode_times=False)
                da = lowpass(da, 13)
                print(run, j, da.std())
                da.plot(c=f'C{i}', ls=ls)
                da.to_netcdf(fn_new)
            
            fn = f'{path_samoc}/SST/TPI1_GMST_dt_raw{tslice}_{run}.nc'
            TPI1 = da = xr.open_dataarray(fn, decode_times=False)
            fn = f'{path_samoc}/SST/TPI2_GMST_dt_raw{tslice}_{run}.nc'
            TPI2 = da = xr.open_dataarray(fn, decode_times=False)
            fn = f'{path_samoc}/SST/TPI3_GMST_dt_raw{tslice}_{run}.nc'
            TPI3 = da = xr.open_dataarray(fn, decode_times=False)
            
            TPI = TPI2 - (TPI1+TPI3)/2
            lowpass(TPI, 13).plot(c='C3', ls=ls)
            lowpass(TPI, 13).to_netcdf(f'{path_samoc}/SST/TPI{tslice}_{run}.nc')
    
    @staticmethod
    def derive_yrly_autocorrelations(run):
        tslices = ['']
        if run=='ctrl':
            tslices.extend(['_100_248', '_151_299'])
            slices = [yrly_ctrl_100_248, yrly_ctrl_151_299]
        elif run=='lpd':
            tslices.extend(['_268_416', '_417_565'])
            slices = [yrly_lpd_268_416, yrly_lpd_417_565]
            
        for j, tslice in enumerate(tslices):
            print(j)
            fn = f'{path_samoc}/SST/SST_GMST_dt_yrly{tslice}_{run}.nc'
            da = xr.open_dataarray(fn, decode_times=False)
#             if j>0:
#                 da = da.sel(time=slice(*slices[j-1]))
            FA = FieldAnalysis(da)
            fn_new = f'{path_samoc}/SST/SST_autocorrelation{tslice}_{run}.nc'
            FA.make_autocorrelation_map(fn_new)
            
    
    @staticmethod
    def make_yrly_regression_files(run, idx):
        """ generate regression files """
        assert idx in ['AMO', 'SOM', 'TPI']
        assert run in ['ctrl', 'lpd', 'had']
        
        tslices = ['']
        if run=='ctrl':
            tslices.extend(['_100_248', '_151_299'])
            slices = [yrly_ctrl_100_248, yrly_ctrl_151_299]
        elif run=='lpd':
            tslices.extend(['_268_416', '_417_565'])
            slices = [yrly_lpd_268_416, yrly_lpd_417_565]
            
        for j, tslice in enumerate(tslices):
            print(j)
            fn = f'{path_samoc}/SST/{idx}{tslice}_{run}.nc'
            index = xr.open_dataarray(fn, decode_times=False)
            fn = f'{path_samoc}/SST/SST_GMST_dt_yrly{tslice}_{run}.nc'
            SST_dt = xr.open_dataarray(fn, decode_times=False)
            fn = f'{path_samoc}/SST/SST_autocorrelation{tslice}_{run}.nc'
            autocorr = xr.open_dataarray(fn, decode_times=False)
            
            xA = xrAnalysis()
            ds = xA.lag_linregress(x=index[7:-7],  # removing filter edge effects
                                   y=SST_dt[7:-7], 
                                   autocorrelation=autocorr,
                                   standardize=True,
                                  )
            ds.to_netcdf(f'{path_samoc}/SST/{idx}_regr{tslice}_{run}.nc')
        print('success')
        return
            

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