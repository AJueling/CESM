import numpy as np
import xarray as xr

from paths import path_samoc
from constants import cp_sw
from timeseries import IterateOutputCESM
from xr_DataArrays import dll_from_arb_da
from xr_DataArrays import generate_xr_DZ, generate_xr_AREA, generate_xr_HTN



def OHC_integrals(domain, run, mask_nr=0):
    """ Ocean Heat Content integration function
    uses either:
    'ocn':      high res. tripolar gridde yrly averaged TEMP/PD fields
    'ocn_rect': low res. yrly averaged ocn_rect TEMP/PD fields
    
    input:
    domain .. (str)
    run    .. (str) 
    
    output:
    ds_new .. xr Dataset
    
    (ocn: )
    (ocn_rect: takes about 3 minutes for all year)
    """
    assert domain in ['ocn', 'ocn_rect']
    assert run in ['ctrl', 'rcp']
    assert type(mask_nr)==int
        
    DZT    = generate_xr_DZ(domain)
    TAREA  = generate_xr_AREA(domain)
    HTN    = generate_xr_HTN(domain)
    
    ds_new = xr.Dataset()
    
    for y,m,file in IterateOutputCESM(domain, run, 'yrly', name='TEMP_PD'):
        print(y)
        t = y*365  # time in days since year 0
        
        ds = xr.open_dataset(file, decode_times=False)
        OHC = (ds.TEMP+2)*ds.PD*cp_sw*1000  # [g/cm^3] to [kg/m^3]
                                            # added 2K to make all temperatures positive
        OHC = OHC.where(OHC<1e9)  # approx max value @ (30 K * 4000 J/K/kg * 1000 kg/m^3)
        OHC = OHC.where(OHC>0)

        ds_g   = assign_time_to_ds(xr_int_global(da=OHC, AREA=TAREA, DZ=DZT),
                                   'OHC_global', t)
        ds_gl  = assign_time_to_ds(xr_int_global_level(da=OHC, AREA=TAREA, DZ=DZT),
                                   'OHC_global_levels', t)
        ds_v   = assign_time_to_ds(xr_int_vertical(da=OHC, DZ=DZT),
                                   'OHC_vertical', t)
        if domain=='ocn':
            ds_z   = assign_time_to_ds(xr_int_zonal_trplr(da=OHC, HTN=HTN, DZ=DZT),
                                       'OHC_zonal', t)
            ds_zl  = assign_time_to_ds(xr_int_zonal_level_trplr(da=OHC, HTN=HTN),
                                       'OHC_zonal_levels', t)
        elif domain=='ocn_rect':
            ds_z   = assign_time_to_ds(xr_int_zonal_rect(da=OHC, HTN=HTN, DZ=DZT),
                                       'OHC_zonal', t)
            ds_zl  = assign_time_to_ds(xr_int_zonal_level_rect(da=OHC, HTN=HTN),
                                       'OHC_zonal_levels', t)

        ds_new = xr.merge([ds_new, ds_g, ds_gl, ds_z, ds_zl, ds_v])
        
        if y==2000 or y==200: break  # for testing only do one year
            
    ds_new.to_netcdf(path=f'{path_samoc}/OHC/OHC_integrals_{run}.nc', mode='w')
    
    return ds_new



def xr_int_global(da, AREA, DZ):
    """ global volume integral *[m^3] """
    (z, lat, lon) = dll_from_arb_da(da)    
    return (da*AREA*DZ).sum(dim=[z, lat, lon])  # 0D


def xr_int_global_level(da, AREA, DZ):
    """ global volume integral *[m^3] """
    (z, lat, lon) = dll_from_arb_da(da)
    return (da*AREA*DZ).sum(dim=[lat, lon])  # 1D (z)


def xr_int_vertical(da, DZ):
    """ vertical integral *[m] """
    (z, lat, lon) = dll_from_arb_da(da)
    return (da*DZ).sum(dim=z)  # 2D (lat, lon)


def xr_int_zonal_rect(da, HTN, DZ):
    """ integral along dpeth and zonal coordinates *[m^2] rectangular grid"""
    (z, lat, lon) = dll_from_arb_da(da)
    return (da*HTN*DZ).sum(dim=[z, lon])  # 1D (lat)


def xr_int_zonal_trplr(da, HTN, DZ):
    """ integral along dpeth and zonal coordinates *[m^2] tripolar grid """
    (z, lat, lon) = dll_from_arb_da(da)
                # proper zonal averaging
    return (da*HTN*DZ).sum(dim=[z, lon])  # 1D (lat)


def xr_int_zonal_level_rect(da, HTN):
    """ zonal integrals for each level *[m] rectangular grid"""
    (z, lat, lon) = dll_from_arb_da(da)
    return (da*HTN).sum(dim=[lon])  # 2D (z, lat)


def xr_int_zonal_level_trplr(da, HTN):
    """ zonal integrals for each level *[m] tripolar grid"""
    (z, lat, lon) = dll_from_arb_da(da)
                # proper zonal averaging
    return (da*HTN).sum(dim=[lon])  # 2D (z, lat)



def assign_time_to_ds(da, name, t):
    """ 
    adds time dimension to xr DataArray, then sets time value to d,
    and then returns xr DataArray as xr dataset
    """
    da = da.expand_dims('time')
    da = da.assign_coords(time=[t])
    ds = da.to_dataset(name=name)
    
    return ds