# from https://gist.github.com/rabernat/bc4c6990eb20942246ce967e6c9c3dbe

import xarray as xr
import numpy as np
from paths import CESM_filename
from timeseries import IterateOutputCESM


def linear_trend(x):
    """ function to compute a linear trend of a timeseries """
    pf = np.polyfit(x.time, x, 1)
    # we need to return a dataarray or else xarray's groupby won't be happy
    return xr.DataArray(pf[0])



def xr_linear_trends_2D(da, dim_names):
    """ calculate linear trend of 2D field in time
    
    input:
    da        .. 3D xr DataArrat with ('lat', 'lon') dimensions
    dim_names .. tuple of 2 strings: lat, lon dimension names
    
    output:
    da_trend  .. slope of linear regression
    """
    (lat, lon) = dim_names
    # stack lat and lon into a single dimension called allpoints
    stacked = da.stack(allpoints=[lat, lon])
    # apply the function over allpoints to calculate the trend at each point
    trend = stacked.groupby('allpoints').apply(linear_trend)
    # unstack back to lat lon coordinates
    da_trend = trend.unstack('allpoints')
    return da_trend



def atm_field_regression(run):
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
    first_year = IterateOutputCESM(domain=domain, run=run, tavg=tavg).year
    first_file = CESM_filename(domain=domain, run=run, y=first_year, m=0, name='T_T850')
    
    da = xr.open_dataset(first_file, decode_times=False)[field][lev,:,:]
    da = da.expand_dims('time')
    da_newer = da.assign_coords(time=[first_year])
    for y, m, file in IterateOutputCESM(domain=domain, run=run, tavg=tavg):
        if y>first_year:
            da_new = xr.open_dataset(file, decode_times=False)[field][lev,:,:]
            da = xr.concat([da, da_new], dim='time')
            print(y)
    
    da_trend = xr_linear_trends_2D(da, ('lat', 'lon'))
    da_trend = da_trend*365*100  # change from [K/day] to [K/century]
    
    return da_trend