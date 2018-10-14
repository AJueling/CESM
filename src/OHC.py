import numpy as np
import xarray as xr
from xr_DataArrays import generate_xr_DZ, generate_xr_AREA, generate_xr_HTN
from paths import path_samoc
from constants import cp_sw

def xr_int_global(da, AREA, DZ):
    return (da*AREA*DZ).sum(dim=['z_t', 'nlat', 'nlon'])  # 0D


def xr_int_global_level(da, AREA, DZ):
    """ """
    return (da*AREA*DZ).sum(dim=['nlat', 'nlon'])  # 1D (z)


def xr_int_zonal(da, HTN, DZ):
    return (da*HTN*DZ).sum(dim=['z_t', 'nlon'])  # 1D (lat)


def xr_int_zonal_level(da, HTN, DZ):
    return (da*HTN*DZ).sum(dim=['nlon'])  # 2D (z, lat)


def xr_int_vertical(da, DZ):
    return (da*DZ).sum(dim='z_t')  # 2D (lat, lon)


def OHC_integrals(run, mask_nr=0):
    
    assert run=='ctrl' or run=='rcp'
    
    DZT = generate_xr_DZ('ocn_lowres_fbc')
    TAREA = generate_xr_AREA('ocn_lowres')
    HTN = generate_xr_HTN('ocn_lowres')
    
    ds_new = xr.Dataset()
    
    # time loop
    ds = xr.open_dataset(f'{path_samoc}/ctrl/ocn_yrly_TEMP_PD_0200.nc', decode_times=False)
    OHC = ds.TEMP*ds.PD*cp_sw
    
    ds_g   = xr_int_global(da=OHC, AREA=TAREA, DZ=DZT).to_dataset(name='OHC_global')
    ds_gl  = xr_int_global_level(da=OHC, AREA=TAREA, DZ=DZT).to_dataset(name='OHC_global_levels')
    ds_z   = xr_int_zonal(da=OHC, HTN=HTN, DZ=DZT).to_dataset(name='OHC_zonal')
    ds_zl  = xr_int_zonal_level(da=OHC, HTN=HTN, DZ=DZT).to_dataset(name='OHC_zonal_levels')
    ds_v   = xr_int_vertical(da=OHC, DZ=DZT).to_dataset(name='OHC_vertical')
    ds_new = xr.merge([ds_new, ds_g, ds_gl, ds_z, ds_zl, ds_v])

    ds_new.to_netcdf(path=f'{path_samoc}/OHC/OHC_integrals_{run}.nc', mode='w')
    
    return ds_new