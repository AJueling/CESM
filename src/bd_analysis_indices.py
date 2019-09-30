import glob
import os
import numpy as np
import xarray as xr
import warnings
import matplotlib.pyplot as plt

from tqdm import tqdm
from eofs.xarray import Eof

from paths import path_samoc, path_prace
from regions import boolean_mask, TPI_masks, mask_box_in_region,\
                    bll_AMO, bll_SOM, bll_TPI1, bll_TPI2, bll_TPI3
from filters import chebychev, lowpass
from timeseries import IterateOutputCESM
from xr_regression import xr_lintrend
from xr_DataArrays import xr_AREA, dll_dims_names

from ab_derivation_SST import DeriveSST, times_ctrl, times_lpd
from ba_analysis_dataarrays import AnalyzeDataArray
from bc_analysis_fields import AnalyzeField

warnings.simplefilter(action='ignore', category=RuntimeWarning)  # to ignore mean of nan message


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
    

    def derive_SST_avg_index(self, run, index, time=None):
        """ generates all area avg indices from detrended SST data """
        assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'had']
        
        if run=='had':
            assert time is None
            domain = 'ocn_had'
            dims = ('latitude', 'longitude')
            ts = ''
            fn_yrly = f'{path_samoc}/SST/SST_GMST_tfdt_yrly_had.nc'
        else:
            assert len(time)==2
            if run in ['lpd', 'lpi']:
                domain = 'ocn_low'
                dims = ('nlat', 'nlon')
            ts = f'_{time[0]}_{time[1]}'
            fn_yrly = f'{path_samoc}/SST/SST_quadratic_pwdt_yrly_{run}{ts}.nc'
        
        if index in ['AMO', 'SOM']:  # yrly data
            SST_yrly = xr.open_dataarray(fn_yrly, decode_times=False)
            if run in ['ctrl', 'rcp']:
                domain = 'ocn'
                dims = ('nlat', 'nlon')
            blats, blons, mask_nr = self.bounding_lats_lons(index)
            MASK = mask_box_in_region(domain=domain, mask_nr=mask_nr, bounding_lats=blats, bounding_lons=blons)
            AREA = xr_AREA(domain=domain).where(MASK)
            SST_index = self.SST_area_average(xa_SST=SST_yrly, AREA=AREA, AREA_index=AREA.sum(), MASK=MASK, dims=dims)
            SST_index.to_netcdf(f'{path_samoc}/SST/{index}_dt_raw_{run}{ts}.nc')
        
        if index=='TPI':  # monthly data
            fn_monthly = f'{path_samoc}/SST/SST_monthly_ds_dt_{run}{ts}.nc'
            SST_monthly = xr.open_dataarray(fn_monthly, decode_times=False)
            if run in ['ctrl', 'rcp']:
                domain = 'ocn_rect'
                dims = ('t_lat', 't_lon')  
            for i, TPI_ in enumerate(['TPI1', 'TPI2', 'TPI3']):
                blats, blons, mask_nr = self.bounding_lats_lons(TPI_)
                MASK = mask_box_in_region(domain=domain, mask_nr=mask_nr, bounding_lats=blats, bounding_lons=blons)
                AREA = xr_AREA(domain=domain).where(MASK)
                SST_index = self.SST_area_average(xa_SST=SST_monthly, AREA=AREA, AREA_index=AREA.sum(), MASK=MASK, dims=dims)
                SST_index.to_netcdf(f'{path_samoc}/SST/{index_}_ds_dt_raw_{run}{ts}.nc')
                if i==0:    TPI = -0.5*SST_index
                elif i==1:  TPI = TPI + SST_index
                elif i==2:  TPI = TPI - 0.5*SST_index
                TPI.to_netcdf(f'{path_samoc}/SST/TPI_ds_dt_raw_{run}{ts}.nc')
            
        return
  
        
    def derive_final_SST_indices(self, run, index, time=None):
        """ processes raw indices: filtering and TPI summation """
        assert run in ['ctrl', 'lpd', 'had']
            
        # AMO & SOM
        if index in ['AMO', 'SOM']:
            fn = f'{path_samoc}/SST/{idx}_{dts}_{dt}_raw_{run}.nc'
            da = xr.open_dataarray(fn, decode_times=False)
            da = DeriveSST().select_time(da, time)
            lowpass(da, 13).to_netcdf(f'{path_samoc}/SST/{idx}_{dts}_{dt}_{run}{ts}.nc')

        elif index in ['TPI']:
            lowpass(TPI, 13*12).to_netcdf(f'{path_samoc}/SST/TPI_{dts}_{dt}_{run}{ts}.nc')
        
        elif index=='PMV':
            for extent in ['38S', 'Eq', '20N']:

        return
    
    
    def derive_yrly_autocorrelations(self, run, tavg, time=None):
        """ autocorrelation maps for detrended SST fields for significance tests """
        assert tavg in ['yrly', 'monthly']
        
        if run=='had':
            if tavg='monthly':  fn = f'{path_samoc}/SST/SST_GMST_tfdt_yrly_had.nc'
            elif tavg='yrly':   fn = f'{path_samoc}/SST/SST_quadratic_pwdt_yrly_had.nc'
            fn_new = f'{path_samoc}/SST/SST_{tavg}_autocorrelation_{run}.nc'
            
        elif run in ['ctrl', 'lpd']:
            if tavg='monthly':  fn = f'{path_samoc}/SST/SST_GMST_tfdt_yrly_{run}_{time[0]}_{time[1]}.nc'
            elif tavg='yrly':   fn = f'{path_samoc}/SST/SST_quadratic_pwdt_yrly_{run}_{time[0]}_{time[1]}.nc'
            fn_new = f'{path_samoc}/SST/SST_{tavg}_autocorrelation_{run}_{time[0]}_{time[1]}.nc'
            
        da = xr.open_dataarray(fn, decode_times=False)
        if time!='full':  da = DeriveSST().select_time(da, time)
        AnalyzeField(da).make_autocorrelation_map(fn_new)
        return
    
    
    def make_yrly_regression_files(self, run, index, time=None):
        """ generate regression files """
        assert run in ['ctrl', 'lpd', 'had']
        
        ts = self.time_slice_string(time)
        dt = self.detrend_string(run)

        fn_acr = f'{path_samoc}/SST/SST_autocorrelation_{run}{ts}.nc'
        autocorr = xr.open_dataarray(fn_acr, decode_times=False)

        if run=='had':
            fn_SST = f'{path_samoc}/SST/SST_GMST_sqdt_yrly_had{ts}.nc'
        elif run in ['ctrl', 'lpd']:
            fn_SST = f'{path_samoc}/SST/SST_quadratic_pwdt_yrly_{run}.nc'
        SST_dt = xr.open_dataarray(fn_SST, decode_times=False)
        
        for idx in ['AMO', 'SOM', 'TPI']:
            fn_idx = f'{path_samoc}/SST/{idx}_{run}.nc'
            index = xr.open_dataarray(fn_idx, decode_times=False)
            if time!='full':  index = DeriveSST().select_time(index, time)

            xA = AnalyzeDataArray()
            ds = xA.lag_linregress(x=index[7:-7],  # removing filter edge effects
                                   y=SST_dt[7:-7], 
                                   autocorrelation=autocorr,
                                   standardize=True,
                                  )
            ds.to_netcdf(f'{path_samoc}/SST/{idx}_regr_{run}{ts}.nc')
        print('success')
        return
    
    
    def EOF_SST_analysis(self, xa, weights, neofs=1, npcs=1, fn=None):
        """ Empirical Orthogonal Function analysis of SST(t,x,y) field; from `SST.py` """
        assert type(xa)==xr.core.dataarray.DataArray
        assert type(weights)==xr.core.dataarray.DataArray
        assert 'time' in xa.dims
        assert np.shape(xa[0,:,:])==np.shape(weights)

        # anomalies by removing time mean
        xa = xa - xa.mean(dim='time')
        # Retrieve the leading EOF, expressed as the covariance between the leading PC
        # time series and the input xa anomalies at each grid point.
        solver = Eof(xa, weights=weights)
        eofs = solver.eofsAsCovariance(neofs=neofs)
        pcs = solver.pcs(npcs=npcs, pcscaling=1)
        if fn!=None:  xr.merge([eofs, pcs]).to_netcdf(fn)
        return eofs, pcs
    
    
    def Pacific_EOF_analysis(self, run, extent, time=None):
        """ """
        # 4:45 for 38S_ctrl, 5:07 for 38S_lpd, 3:42 for 38S_had : total 11:08
        assert run in ['ctrl', 'lpd', 'had']
        assert extent in ['38S', 'Eq', '20N']

        if run=='ctrl':
            monthly_fns  = [f'{path_prace}/SST/SST_monthly_ds_dt_{extent}_ctrl_{time[0]}_{time[1]}.nc' for time in times_ctrl]
            EOF_fns      = [f'{path_prace}/SST/PMV_EOF_{extent}_ctrl_{time[0]}_{time[1]}.nc' for time in times_ctrl]
            yrly_times   = xr.open_dataarray(f'{path_prace}/SST/SST_yrly_rect_ctrl.nc', decode_times=False).time
            domain       = 'ocn_rect'
        elif run=='lpd':
            monthly_fns  = [f'{path_prace}/SST/SST_monthly_ds_dt_{extent}_lpd_{time[0]}_{time[1]}.nc' for time in times_lpd]
            EOF_fns      = [f'{path_prace}/SST/PMV_EOF_{extent}_lpd_{time[0]}_{time[1]}.nc' for time in times_lpd]
            yrly_times   = xr.open_dataarray(f'{path_prace}/SST/SST_yrly_lpd.nc', decode_times=False).time
            domain       = 'ocn_low'
        elif run=='had':
            assert time is None
            monthly_fns  = [f'{path_prace}/SST/SST_monthly_ds_dt_{extent}_had.nc']
            EOF_fns      = [f'{path_prace}/SST/PMV_EOF_{extent}_had.nc']
            yrly_times   = xr.open_dataarray(f'{path_prace}/SST/SST_yrly_had.nc', decode_times=False).time
            domain       = 'ocn_had'
            
        area = xr.open_dataarray(f'{path_prace}/geometry/AREA_{extent}_{domain}.nc')
        (d, lat, lon) = dll_dims_names(domain)
            
        for k, fn in tqdm(enumerate(monthly_fns)):
            if k<=2:  continue
                
            fn = monthly_fns[k]
            fn_EOF = EOF_fns[k]
            da = xr.open_dataarray(fn, decode_times=False)
            eofs, pcs = self.EOF_SST_analysis(xa=da, weights=area, neofs=1, npcs=1, fn=None)
            cov = eofs.mean(dim=[lat,lon])
            if cov<0:  factor=-1
            else:      factor= 1
            if run in ['lpd', 'had'] and extent=='20N':  factor = factor*-1
            
            ds = xr.merge([eofs*factor, pcs*factor])
            ds.to_netcdf(fn_EOF)
            
        return ds, yrly
    
        
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
            
    
    def time_slice_string(self, time):
        """ string for time subset """
        if time=='full':         ts = ''
        elif type(time)==tuple:  ts = f'_{time[0]}_{time[1]}'
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