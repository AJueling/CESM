import numpy as np
import xarray as xr

from paths import path_samoc
from regions import boolean_mask, regions_dict
from constants import cp_sw
from timeseries import IterateOutputCESM
from xr_integrate import xr_int_global, xr_int_global_level, xr_int_vertical,\
                         xr_int_zonal, xr_int_zonal_level
from xr_DataArrays import xr_DZ, xr_AREA, xr_HTN, xr_LATS, dll_from_arb_da


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
    DZT  = xr_DZ(domain).where(MASK)
    AREA = xr_AREA(domain).where(MASK)
    HTN  = xr_HTN(domain).where(MASK)
    LATS = xr_LATS(domain).where(MASK)
    
    ds_new = xr.Dataset()
    
    for y,m,file in IterateOutputCESM(domain, run, 'yrly', name='TEMP_PD'):
        print(y, file)
        t   = y*365  # time in days since year 0
        ds  = xr.open_dataset(file, decode_times=False)
        OHC = ds.TEMP*ds.PD*cp_sw*1000  # [g/cm^3] to [kg/m^3]
        OHC = OHC.where(MASK)

        # xr DataArrays
        da_g  = xr_int_global(da=OHC, AREA=AREA, DZ=DZT)
        da_gl = xr_int_global_level(da=OHC, AREA=AREA, DZ=DZT)
        da_v  = xr_int_vertical(da=OHC, DZ=DZT)
        da_z  = xr_int_zonal(da=OHC, HTN=HTN, LATS=LATS, AREA=AREA, DZ=DZT)
        da_zl = xr_int_zonal_level(da=OHC, HTN=HTN, LATS=LATS, AREA=AREA, DZ=DZT)

        # xr Datasets
        ds_g  = t2ds(da_g , 'OHC_global'       , t)
        ds_gl = t2ds(da_gl, 'OHC_global_levels', t)
        ds_v  = t2ds(da_v , 'OHC_vertical'     , t)
        ds_z  = t2ds(da_z , 'OHC_zonal'        , t)
        ds_zl = t2ds(da_zl, 'OHC_zonal_levels' , t)
        
        ds_new = xr.merge([ds_new, ds_g, ds_gl, ds_z, ds_zl, ds_v])
        
        if y==2000 or y==200: break  # for testing only
            
    path_out = f'{path_samoc}/OHC/OHC_integrals_{regions_dict[mask_nr]}_{run}.nc'
    print(f'output: {path_out}')
    ds_new.to_netcdf(path=path_out, mode='w')
    
    return ds_new



def t2ds(da, name, t):
    """ 
    adds time dimension to xr DataArray, then sets time value to t,
    and then returns xr dataset
    """
    da = da.expand_dims('time')
    da = da.assign_coords(time=[t])
    ds = da.to_dataset(name=name)
    
    return ds