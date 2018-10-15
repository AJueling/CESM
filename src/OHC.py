import numpy as np
import xarray as xr
from xr_DataArrays import generate_xr_DZ, generate_xr_AREA, generate_xr_HTN
from paths import path_samoc
from constants import cp_sw
from timeseries import IterateOutputCESM

def xr_int_global(da, AREA, DZ):
    return (da*AREA*DZ).sum(dim=['depth_t', 'lat', 'lon'])  # 0D


def xr_int_global_level(da, AREA, DZ):
    """ global volume integral *[m^3] """
    return (da*AREA*DZ).sum(dim=['lat', 'lon'])  # 1D (z)


def xr_int_zonal(da, HTN, DZ):
    """ integral along dpeth and zonal coordinates *[m^2] """
    return (da*HTN*DZ).sum(dim=['depth_t', 'lon'])  # 1D (lat)


def xr_int_zonal_level(da, HTN):
    """ zonal integrals for each level *[m] """
    return (da*HTN).sum(dim=['lon'])  # 2D (z, lat)


def xr_int_vertical(da, DZ):
    """ vertical integral *[m] """
    return (da*DZ).sum(dim='depth_t')  # 2D (lat, lon)


def OHC_integrals(run, mask_nr=0):
    """ Ocean Heat Content integration function
    uses low resolution 'rect' yrly averaged ocn TEMP/PD fields
    
    input:
    run .. (str) 
    
    output:
    ds_new .. xr Dataset
    
    (takes about 3 minutes)
    """
    assert run=='ctrl' or run=='rcp'
    
    DZT   = generate_xr_DZ('ocn_rect')
    TAREA = generate_xr_AREA('ocn_rect')
    HTN   = generate_xr_HTN('ocn_rect')
    
    ds_new = xr.Dataset()
    
    # time loop
    for y,m,file in IterateOutputCESM('ocn_rect', run, 'yrly', name='TEMP_PD'):
        print(y)
        d = y*365  # time in days since year 0
        ds = xr.open_dataset(file, decode_times=False)
        OHC = (ds.TEMP+2)*ds.PD*cp_sw*1000  # [g/cm^3] to [kg/m^3]
                                            # added 2K to make is all positive
        OHC = OHC.where(OHC<1e9)  # approx max value @ (30 K * 4000 J/K/kg * 1000 kg/m^3)
        OHC = OHC.where(OHC>0)
#         OHC = OHC.rename({'t_lat': 'lat', 't_lon': 'lon'})

        ds_g   = assign_time_to_ds(xr_int_global(da=OHC, AREA=TAREA, DZ=DZT),
                                   'OHC_global', d)
        ds_gl  = assign_time_to_ds(xr_int_global_level(da=OHC, AREA=TAREA, DZ=DZT),
                                   'OHC_global_levels', d)
        ds_z   = assign_time_to_ds(xr_int_zonal(da=OHC, HTN=HTN, DZ=DZT),
                                   'OHC_zonal', d)
        ds_zl  = assign_time_to_ds(xr_int_zonal_level(da=OHC, HTN=HTN),
                                   'OHC_zonal_levels', d)
        ds_v   = assign_time_to_ds(xr_int_vertical(da=OHC, DZ=DZT),
                                   'OHC_vertical', d)

        ds_new = xr.merge([ds_new, ds_g, ds_gl, ds_z, ds_zl, ds_v])
        
    ds_new.to_netcdf(path=f'{path_samoc}/OHC/OHC_integrals_{run}.nc', mode='w')
    
    return ds_new


def assign_time_to_ds(da, name, d):
    """ 
    adds time dimension to xr DataArray, then sets time value to d,
    and then returns xr DataArray as xr dataset
    """
    da = da.expand_dims('time')
    da = da.assign_coords(time=[d])
    ds = da.to_dataset(name=name)
    
    return ds