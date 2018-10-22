import os
import numpy as np
import xarray as xr

from paths import path_samoc
from regions import boolean_mask, regions_dict
from constants import cp_sw, rho_sw
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
    DZT  = xr_DZ(domain)
    AREA = xr_AREA(domain)
    HTN  = xr_HTN(domain)
    LATS = xr_LATS(domain)
    
    for k in range(42):
        DZT[k,:,:]  = DZT[k,:,:].where(MASK)
    AREA = AREA.where(MASK)
    HTN  = HTN.where(MASK)
    LATS = LATS.where(MASK)
    
#     DZT.TLONG  = DZT.TLONG.round(decimals=2)
#     AREA.TLONG = AREA.TLONG.round(decimals=2)
    HTN['TLAT']   = HTN['TLAT'].round(decimals=2)
    HTN['TLONG']  = HTN['TLONG'].round(decimals=2)
    LATS['TLAT']  = LATS['TLAT'].round(decimals=2)
    LATS['TLONG'] = LATS['TLONG'].round(decimals=2)
    
    file_out = f'{path_samoc}/OHC/OHC_integrals_{run}_{regions_dict[mask_nr]}.nc'
    if os.path.isfile(file_out):  os.remove(file_out)
    print(f'output: {file_out}')
    first_yr = 0
    
    for y,m,file in IterateOutputCESM(domain, run, 'yrly', name='TEMP_PD'):
        if first_yr==0: 
            first_yr=y
            print('first year: ', first_yr)
        print(y, file)
#         if y in [2001,2002,2003]: continue
        t   = y*365  # time in days since year 0
        ds  = xr.open_dataset(file, decode_times=False)
        ds['TLAT'] = ds['TLAT'].round(decimals=2)
        ds['TLONG'] = ds['TLONG'].round(decimals=2)
        
        if ds.PD[0,1000,1000].round(decimals=0)==0:
            print(ds.PD[0,1000,1000])
            ds['PD'] = ds['PD']+rho_sw
        elif ds.PD[0,1000,1000].round(decimals=0)==1:
            print(ds.PD[0,1000,1000])
            pass
        else: 'density is neither close to 0 or 1'
            
        
        OHC = ds.TEMP*ds.PD*cp_sw*1000  # [g/cm^3] to [kg/m^3]
        OHC = OHC.where(MASK)

        # xr DataArrays
        da_g  = xr_int_global(da=OHC, AREA=AREA, DZ=DZT)
        da_gl = xr_int_global_level(da=OHC, AREA=AREA, DZ=DZT)
        da_v  = xr_int_vertical(da=OHC, DZ=DZT)
        da_z  = xr_int_zonal(da=OHC, HTN=HTN, LATS=LATS, AREA=AREA, DZ=DZT)
        da_zl = xr_int_zonal_level(da=OHC, HTN=HTN, LATS=LATS, AREA=AREA, DZ=DZT)

        ds.close()
        
        # xr Datasets
        ds_g  = t2ds(da_g , 'OHC_global'       , t)
        ds_gl = t2ds(da_gl, 'OHC_global_levels', t)
        ds_v  = t2ds(da_v , 'OHC_vertical'     , t)
        ds_z  = t2ds(da_z , 'OHC_zonal'        , t)
        ds_zl = t2ds(da_zl, 'OHC_zonal_levels' , t)
        
        if y==first_yr: 
            ds_new = xr.merge([ds_g, ds_gl, ds_z, ds_zl, ds_v])
        elif y>first_yr:
            ds_temp = xr.open_dataset(file_out, decode_times=False)
            ds_new = xr.merge([ds_temp, ds_g, ds_gl, ds_z, ds_zl, ds_v])
            ds_temp.close()
        ds_new.to_netcdf(path=file_out, mode='w')
        ds_new.close()
        
#         if y==2001 or y==201: break  # for testing only
    
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