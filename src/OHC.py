import os
import numpy as np
import xarray as xr

from grid import create_dz_mean
from paths import path_samoc
from regions import boolean_mask, regions_dict
from constants import cp_sw, rho_sw, km
from timeseries import IterateOutputCESM
from xr_integrate import xr_int_global, xr_int_global_level, xr_int_vertical,\
                         xr_int_zonal, xr_int_zonal_level
from xr_DataArrays import xr_DZ, xr_AREA, xr_HTN, xr_LATS, dll_from_arb_da
from xr_regression import xr_linear_trend


def OHC_integrals(domain, run, mask_nr=0):
    """ Ocean Heat Content integration function
    uses either:
    'ocn':      high res. tripolar gridded yrly averaged TEMP/PD fields
    'ocn_rect': low res. yrly averaged ocn_rect TEMP/PD fields
    
    input:
    domain .. (str)
    run    .. (str) 
    
    output:
    ds_new .. xr Dataset
    
    (ocn:      takes about 45 seconds per year: 70 yrs approx 55 mins)
    (ocn_rect: takes about  3 seconds per year: 70 yrs approx 3 mins)
    """
    assert domain in ['ocn', 'ocn_rect']
    assert run in ['ctrl', 'rcp']
    assert type(mask_nr)==int
    
    MASK = boolean_mask(domain, mask_nr)
    DZT  = xr_DZ(domain)
    AREA = xr_AREA(domain)
    HTN  = xr_HTN(domain)
    LATS = xr_LATS(domain)
    
    for k in range(42):
        DZT[k,:,:]  = DZT[k,:,:].where(MASK)
    AREA = AREA.where(MASK)
    HTN  = HTN.where(MASK)
    LATS = LATS.where(MASK)
    
    HTN['TLAT']   = HTN['TLAT'].round(decimals=2)
    HTN['TLONG']  = HTN['TLONG'].round(decimals=2)
    LATS['TLAT']  = LATS['TLAT'].round(decimals=2)
    LATS['TLONG'] = LATS['TLONG'].round(decimals=2)
    
    for y,m,file in IterateOutputCESM(domain, run, 'yrly', name='TEMP_PD'):
        print(y, file)
        t   = y*365  # time in days since year 0, for consistency with CESM date output
        ds  = xr.open_dataset(file, decode_times=False)
        ds['TLAT'] = ds['TLAT'].round(decimals=2)
        ds['TLONG'] = ds['TLONG'].round(decimals=2)
        
        if ds.PD[0,1000,1000].round(decimals=0)==0:
            ds['PD'] = ds['PD']+rho_sw
        elif ds.PD[0,1000,1000].round(decimals=0)==1:
            pass
        else: 'density is neither close to 0 or 1'
            
        
        OHC = ds.TEMP*ds.PD*cp_sw*1000  # [g/cm^3] to [kg/m^3]
        OHC = OHC.where(MASK)
        
        OHC_DZT = OHC*DZT

        # xr DataArrays
        da_g  = xr_int_global(da=OHC, AREA=AREA, DZ=DZT)
        da_gl = xr_int_global_level(da=OHC, AREA=AREA, DZ=DZT)
        da_v  = OHC_DZT.sum(dim='z_t') #xr_int_vertical(da=OHC, DZ=DZT)
        da_va = OHC_DZT.isel(z_t=slice(0, 9)).sum(dim='z_t')  # above 100 m
        da_vb = OHC_DZT.isel(z_t=slice(9,42)).sum(dim='z_t')  # below 100 m
        da_z  = xr_int_zonal(da=OHC, HTN=HTN, LATS=LATS, AREA=AREA, DZ=DZT)
        da_zl = xr_int_zonal_level(da=OHC, HTN=HTN, LATS=LATS, AREA=AREA, DZ=DZT)

        ds.close()
        
        # xr Datasets
        ds_g  = t2ds(da_g , 'OHC_global'             , t)
        ds_gl = t2ds(da_gl, 'OHC_global_levels'      , t)
        ds_v  = t2ds(da_v , 'OHC_vertical'           , t)
        ds_va = t2ds(da_va, 'OHC_vertical_above_100m', t)
        ds_vb = t2ds(da_vb, 'OHC_vertical_below_100m', t)
        ds_z  = t2ds(da_z , 'OHC_zonal'              , t)
        ds_zl = t2ds(da_zl, 'OHC_zonal_levels'       , t)
        
        
        file_out = f'{path_samoc}/OHC/OHC_integrals_{regions_dict[mask_nr]}_{run}_{y}.nc'
        print(f'output: {file_out}')
        
        ds_new = xr.merge([ds_g, ds_gl, ds_z, ds_zl, ds_v, ds_va, ds_vb])
        ds_new.to_netcdf(path=file_out, mode='w')
        ds_new.close()
        
        if y==2002 or y==102: break  # for testing only
    
    # combining yearly files
    file_out = f'{path_samoc}/OHC/OHC_integrals_{regions_dict[mask_nr]}_{run}.nc'
#    if os.path.isfile(file_out):  os.remove(file_out)
    combined = xr.open_mfdataset(f'{path_samoc}/OHC/OHC_integrals_{regions_dict[mask_nr]}_{run}_*.nc',
                                 concat_dim='time',
                                 autoclose=True,
                                 coords='minimal')
    combined.to_netcdf(f'{path_samoc}/OHC/OHC_integrals_{regions_dict[mask_nr]}_{run}.nc')
    
    return



def t2da(da, t):
    """adds time dimension to xr DataArray, then sets time value to t"""
    da = da.expand_dims('time')
    da = da.assign_coords(time=[t])
    return da


def t2ds(da, name, t):
    """ 
    adds time dimension to xr DataArray, then sets time value to t,
    and then returns as array in xr dataset
    """
    da = t2da(da, t)
    ds = da.to_dataset(name=name)
    
    return ds



def trend_global_levels(ds):
    """ trend of OHC per level as [J/m/y] """
    dz_mean = create_dz_mean(domain='ocn')
    
    return xr_linear_trend(ds.OHC_global_levels/dz_mean)*365



def OHC_detrend_levels(da, detrend='lin'):
    """ """
    assert detrend in ['lin', 'quad']
    assert 'time' in da.coords
    
    dz_mean = create_dz_mean(domain='ocn')
    n = len(da.time)
    times = np.arange(n)
    levels_trend = np.zeros((n, km))
    
    if detrend=='lin':
        lfit_par  = np.zeros((2, km))
        for k in range(km):
            lfit_par[:,k] = np.polyfit(times, da[:,k]/dz_mean[k], 1)

        lin_fit  = np.zeros((n, km))
        for t in range(n):
            lin_fit[t,:]  = lfit_par[0,:]*t + lfit_par[1,:]

        levels_trend = ((da[:,:]/dz_mean - lin_fit ))

    elif detrend=='quad':
        qfit_par = np.zeros((3, km))
        for k in range(km):
            qfit_par[:,k] = np.polyfit(times, da[:,k]/dz_mean[k], 2)
        
        quad_fit = np.zeros((n, km))
        for t in range(n):
            quad_fit[t,:] = qfit_par[0,:]*t**2 + qfit_par[1,:]*t + qfit_par[2,:]
        
        levels_trend = ((da[:,:]/dz_mean - quad_fit))

    return levels_trend



def OHC_vert_diff_mean_rm(ds, run):

    assert run in ['ctrl', 'rcp']
    
    for suffix in ['']:#, '_above_100m','_below_100m']:
        assert f'OHC_vertical{suffix}' in ds
    
        OHC_vert_diff = ds[f'OHC_vertical{suffix}']-ds[f'OHC_vertical{suffix}'].shift(time=1)

        OHC_vert_diff_mean = OHC_vert_diff.mean(dim='time')  # 1 min
        OHC_vert_diff_rm   = OHC_vert_diff.rolling({'time':10}, center=True).mean(dim='time')

        OHC_vert_diff_mean.to_netcdf(f'{path_samoc}/OHC/OHC_vert{suffix}_diff_mean_{run}.nc' )
        OHC_vert_diff_rm  .to_netcdf(f'{path_samoc}/OHC/OHC_vert{suffix}_diff_rm_{run}.nc'  )

    return