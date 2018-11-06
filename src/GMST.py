import os
import numpy as np
import xarray as xr

from paths import path_results
from constants import abs_zero
from timeseries import IterateOutputCESM
from xr_DataArrays import xr_AREA



def GMST_timeseries(run):
    """ builds a timesries of the GMST and saves it to a netCDF
    
    input:
    run    .. (str) ctrl or cp
    
    output:
    ds_new .. xr Dataset containing GMST and T_zonal
    """
    
    domain = 'atm'
    tavg   = 'yrly'
    name   = 'T_T850_U_V'
    
    iterator = IterateOutputCESM(domain=domain, run=run, tavg=tavg, name=name)
    ny       = len(iterator)
    iterator = IterateOutputCESM(domain=domain, run=run, tavg=tavg, name=name)
    first_yr = iterator.year
    years    = np.arange(ny) + first_yr
        
    AREA       = xr_AREA('atm')
    AREA_lat   = AREA.sum(dim='lon')
    AREA_total = AREA.sum(dim=('lat','lon'))
    
    iterator = IterateOutputCESM(domain=domain, run=run, tavg=tavg, name=name)
    for i, (y, m, file) in enumerate(iterator):
        print(y)
        assert os.path.exists(file)
        ds = xr.open_dataset(file, decode_times=False)
        
        if i==0:  # create new xr Dataset
            lats = ds.lat.values
            ds_new = xr.Dataset()
            ds_new['GMST']    = xr.DataArray(data=np.zeros((ny)),
                                             coords={'year': years},
                                             dims=('year'))
            ds_new['T_zonal'] = xr.DataArray(data=np.zeros((ny, len(lats))), 
                                             coords={'year': years, 'lat': lats},
                                             dims=('year', 'lat'))
            
        ds_new['GMST'][i]      = (ds['T'][-1,:,:]*AREA).sum(dim=('lat','lon'))/AREA_total
        ds_new['T_zonal'][i,:] = (ds['T'][-1,:,:]*AREA).sum(dim='lon')/AREA_lat
        
    for field in ['GMST', 'T_zonal']:  # [K] to [degC]
        ds_new[field] = ds_new[field] + abs_zero
        
    ds_new.to_netcdf(path=f'{path_results}/GMST/GMST_{run}.nc', mode='w')
    
    return ds_new



def GMST_regression(run):
    """ calculates the surface temperature trends
    
    input:
    run      .. (str) ctrl or rcp
    
    output:
    da_trend .. 2D xr DataArray with linear trends
    
    (takes about 4:30 minutes)
    """
    
    field = 'T'
    lev = -1
    
    domain = 'atm'
    tavg   = 'yrly'
    name   = 'T_T850_U_V'
    first_year = IterateOutputCESM(domain=domain, run=run, tavg=tavg, name=name).year
    first_file = CESM_filename(domain=domain, run=run, y=first_year, m=0, name=name)
    
    da = xr.open_dataset(first_file, decode_times=False)[field][lev,:,:]
    da = da.expand_dims('time')
    for y, m, file in IterateOutputCESM(domain=domain, run=run, tavg=tavg, name=name):
        if y>first_year:
            da_new = xr.open_dataset(file, decode_times=False)[field][lev,:,:]
            da = xr.concat([da, da_new], dim='time')
            print(y)
    
    da_trend = xr_linear_trends_2D(da, ('lat', 'lon'))
    da_trend = da_trend*365*100  # change from [K/day] to [K/century]
    
    return da_trend