import glob
import os
import xarray as xr
import matplotlib.pyplot as plt

from SST import SST_index
from paths import path_samoc
from timeseries import IterateOutputCESM
from xr_regression import xr_lintrend

class GenerateSSTFields(object):
    """"""
    def __init__(self):
        return
    
    @staticmethod
    def remove_superfluous_files(fn):
        for x in glob.glob(fn):
            os.remove(x) 
    
    @staticmethod
    def generate_yrly_SST_files(run):
        """ generate the SST data files from TEMP_PD yaryl averaged files """
        # ca. 4:30 min for ctrl/rcp, 1:25 for lpi
        # stacking files into one xr DataArray object
        for i, (y,m,s) in enumerate(IterateOutputCESM('ocn', run, 'yrly', name='TEMP_PD')):
            print(y)
            da = xr.open_dataset(s, decode_times=False).TEMP[0,:,:]
            da = da.drop(['z_t', 'ULONG', 'ULAT'])
            da['TLAT' ] = da['TLAT' ].round(decimals=2)
            da['TLONG'] = da['TLONG'].round(decimals=2)
            del da.encoding["contiguous"]
            ds = t2ds(da=da, name='SST', t=int(round(da.time.item())))
            ds.to_netcdf(path=f'{path_samoc}/SST/SST_yrly_{run}_{y}.nc', mode='w')

        combined = xr.open_mfdataset(f'{path_samoc}/SST/SST_yrly_{run}_*.nc',
                                     concat_dim='time', autoclose=True, coords='minimal')
        combined.to_netcdf(f'{path_samoc}/SST/SST_yrly_{run}.nc')

        GenerateSSTFields.remove_superfluous_files(f'{path_samoc}/SST/SST_yrly_{run}_*.nc')
            
            
    @staticmethod
    def generate_monthly_SST_files(run):
        """ concatonate monthly files, ocn_rect for high res runs"""
        # 8 mins for 200 years of ctrl
        if run in ['ctrl', 'rcp']:  domain = 'ocn_rect'
        elif run in ['lpd', 'lpi']:  domain = 'ocn_low'
            
        for y,m,s in IterateOutputCESM(domain=domain, tavg='monthly', run=run):
            if run in ['ctrl', 'rcp']:
                xa = xr.open_dataset(s, decode_times=False).TEMP[0,:,:]
            if run in ['lpd', 'lpi']:
                xa = xr.open_dataset(s, decode_times=False).TEMP[0,0,:,:]
            if m==1:
                print(y)
                xa_out = xa.copy()    
            else:
                xa_out = xr.concat([xa_out, xa], dim='time')
            if m==12:
                xa_out.to_netcdf(f'{path_samoc}/SST/SST_monthly_{run}_{y}.nc')
                        
        combined = xr.open_mfdataset(f'{path_samoc}/SST/SST_monthly_{run}_*.nc',
                                     concat_dim='time', decode_times=False)
        combined.to_netcdf(f'{path_samoc}/SST/SST_monthly_{run}.nc')
        combined.close()

        GenerateSSTFields.remove_superfluous_files(f'{path_samoc}/SST/SST_monthly_{run}_*.nc')
            
            
    @staticmethod
    def generate_monthly_regional_SST_files(run):
        """"""
            # 6 mins for 200 years of ctrl
        for i, r in enumerate(['Pac_38S', 'Pac_Eq', 'Pac_20N']):
            latS = [-38, 0, 20][i]
            lonE = [300, 285, 255][i]
            
            if run  in ['ctrl', 'rcp']:  domain= 'ocn_rect'
            elif run in ['lpd', 'lpi']:  domain = 'ocn_low'
                
            Pac_MASK = mask_box_in_region(domain=domain, mask_nr=2,
                                          bounding_lats=(latS,68),
                                          bounding_lons=(110,lonE))
            if run in ['ctrl', 'rcp']:
                Pac_MASK = Pac_MASK.where(Pac_MASK.t_lon+1/.6*Pac_MASK.t_lat<333,0)
            NPac_area = xr_AREA(domain).where(Pac_MASK, drop=True)
            
            for y,m,s in IterateOutputCESM(domain=domain, tavg='monthly', run=run):
                xa = xr.open_dataset(s, decode_times=False).TEMP[0,:,:].where(Pac_MASK, drop=True)
                if m==1:
                    print(y)
                    xa_out = xa.copy()    
                else:
                    xa_out = xr.concat([xa_out, xa], dim='time')
                if m==12:
                    xa_out.to_netcdf(f'{path_samoc}/SST/SST_monthly_{r}_{run}_{y}.nc')
                            
            combined = xr.open_mfdataset(f'{path_samoc}/SST/SST_monthly_{r}_{run}_*.nc',
                                         concat_dim='time', decode_times=False)
            combined.to_netcdf(f'{path_samoc}/SST/SST_monthly_{r}_{run}.nc')
            combined.close()
            
            GenerateSSTFields.remove_superfluous_files(f'{path_samoc}/SST/SST_monthly_{r}_{run}_*.nc')
            # # remove yearly files
            
    @staticmethod
    def generate_monthly_mock_linear_GMST_files(run):
        """ generates xr DataArray with same time coordinates as monthly SST fields
        which contains the linear fit to the mean SST as a stand in for the missing
        monthly GMST
        """
        assert run in ['ctrl', 'lpd']
        if run=='ctrl':  dims = ['t_lon', 't_lat']
        elif run=='lpd':  dims = ['nlon', 'nlat']
        da = xr.open_dataarray(f'{path_samoc}/SST/SST_monthly_{run}.nc', decode_times=False)
        da_new = xr_lintrend(da.mean(dim=dims, skipna=True, keep_attrs=True))
        da_new.name = 'GMST'
        da_new.attrs = {'Note':'This is the linear trend of the SST evolution, not GMST'}
        da_new.to_dataset().to_netcdf(f'{path_samoc}/GMST/GMST_monthly_{run}.nc')
        

yrly_ctrl_100_248 = (36531,90551)
yrly_ctrl_151_299 = (55146,109166)
yrly_lpd_268_416  = (97851,151871)
yrly_lpd_417_565  = (152236,206256)

from filters import chebychev, lowpass
from analysis import FieldAnalysis, xrAnalysis

class DeriveSSTIndices(object):
    def __init__():
        return
    
    @staticmethod
    def derive_all_SST_avg_indices(run):
        """"""
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