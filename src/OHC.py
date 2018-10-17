import numpy as np
import xarray as xr

from paths import path_samoc
from regions import boolean_mask
from constants import cp_sw
from timeseries import IterateOutputCESM
from xr_integrate import xr_int_global, xr_int_global_level, xr_int_vertical,\
                         xr_int_zonal, xr_int_zonal_level
from xr_DataArrays import generate_xr_DZ, generate_xr_AREA, generate_xr_HTN,\
                          dll_from_arb_da



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
        
    DZT  = generate_xr_DZ(domain)
    AREA = generate_xr_AREA(domain)
#     HTN  = generate_xr_HTN(domain)
    MASK = boolean_mask(domain, mask_nr)
    
    ds_new = xr.Dataset()
    
    for y,m,file in IterateOutputCESM(domain, run, 'yrly', name='TEMP_PD'):
        print(y)
        t = y*365  # time in days since year 0
        
        ds = xr.open_dataset(file, decode_times=False)
        OHC = ds.TEMP*ds.PD*cp_sw*1000  # [g/cm^3] to [kg/m^3]
#                                             added 2K to make all temperatures positive
        OHC = OHC.where(MASK)
        OHC = OHC.where(OHC<1e9)  # approx max value @ (30 K * 4000 J/K/kg * 1000 kg/m^3)
#         OHC = OHC.where(OHC>0)

#         ds_g  = t2da(xr_int_global(da=OHC, AREA=AREA, DZ=DZT)      , 'OHC_global'       , t)
#         ds_gl = t2da(xr_int_global_level(da=OHC, AREA=AREA, DZ=DZT), 'OHC_global_levels', t)
#         ds_v  = t2da(xr_int_vertical(da=OHC, DZ=DZT)               , 'OHC_vertical'     , t)
        ds_z  = t2da(xr_int_zonal(da=OHC, AREA=AREA, DZ=DZT)       , 'OHC_zonal'        , t)
        ds_zl = t2da(xr_int_zonal_level(da=OHC, AREA=AREA)         , 'OHC_zonal_levels' , t)

        ds_new = xr.merge([ds_new, ds_z, ds_zl])
#         ds_new = xr.merge([ds_new, ds_g, ds_gl, ds_z, ds_zl, ds_v])
        
        if y==2000 or y==201: break  # for testing only
            
    ds_new.to_netcdf(path=f'{path_samoc}/OHC/OHC_integrals_{run}.nc', mode='w')
    
    return ds_new


def t2da(da, name, t):
    """ 
    adds time dimension to xr DataArray, then sets time value to d,
    and then returns xr DataArray as xr dataset
    """
    da = da.expand_dims('time')
    da = da.assign_coords(time=[t])
    ds = da.to_dataset(name=name)
    
    return ds