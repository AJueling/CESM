# from https://gist.github.com/rabernat/bc4c6990eb20942246ce967e6c9c3dbe

import xarray as xr
import numpy as np
from paths import CESM_filename
from timeseries import IterateOutputCESM


def xr_linear_trend(x):
    """ function to compute a linear trend of a timeseries """
#     print(type(x))  # DataArray
#     varx = x.time.values
#     vary = x.values
#     idx = np.isfinite(varx) & np.isfinite(vary)
#     mask = ~np.isnan(varx) & ~np.isnan(vary)
    pf = np.polynomial.polynomial.polyfit(x.time, x, 1)
    # we need to return a dataarray or else xarray's groupby won't be happy
    return xr.DataArray(pf[1])



def xr_linear_trends_2D(da, dim_names):
    """ calculate linear trend of 2D field in time
    
    input:
    da        .. 3D xr DataArrat with ('lat', 'lon') dimensions
    dim_names .. tuple of 2 strings: e.g. lat, lon dimension names
    
    output:
    da_trend  .. slope of linear regression
    """
    (dim1, dim2) = dim_names
    # stack lat and lon into a single dimension called allpoints
    stacked = da.stack(allpoints=[dim1, dim2])
    # apply the function over allpoints to calculate the trend at each point
    trend = stacked.groupby('allpoints').apply(xr_linear_trend)
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
#     da_newer = da.assign_coords(time=[first_year])
    for y, m, file in IterateOutputCESM(domain=domain, run=run, tavg=tavg):
        if y>first_year:
            da_new = xr.open_dataset(file, decode_times=False)[field][lev,:,:]
            da = xr.concat([da, da_new], dim='time')
            print(y)
    
    da_trend = xr_linear_trends_2D(da, ('lat', 'lon'))
    da_trend = da_trend*365*100  # change from [K/day] to [K/century]
    
    return da_trend



def zonal_trend(ds):
    """ calculated tends for lat_bin field
    
    input:
    ds    .. xr Dataset of OHC fields
             ds.OHC_zonal [J m^-1]
    
    output:
    trend .. 1D xr DataArray in [J m^-1 year^-1]
    """    
    assert 'OHC_zonal' in ds
    assert 'time' in ds
    (min_lat, max_lat) = find_valid_domain(ds.OHC_zonal)
    trend = xr_linear_trend(ds.OHC_zonal[:,min_lat:max_lat+1])*365
    trend = trend.rename({'dim_0': 'TLAT_bins'})
    trend = trend.assign_coords(TLAT_bins=ds.TLAT_bins[min_lat:max_lat+1])
    return trend



def zonal_levels_trend(ds):
    """ calculated trend for depth-lat_bin field
    
    input:
    ds .. xr Dataset
          ds.OHC_zonal_levels [J m^-2]
    
    output:
    trend .. xr DataArray containing trends of OHC [J m^-2 year^-1]
    """
    
    assert 'OHC_zonal_levels' in ds
    
    dn = ('TLAT_bins', 'z_t')
    (min_lat, max_lat, max_depth) = find_valid_domain(ds.OHC_zonal_levels)
    print(min_lat, max_lat)
    trend = xr_linear_trends_2D(ds.OHC_zonal_levels[:,:max_depth+1,min_lat:max_lat+1], dim_names=dn)*365  #  to [yr^-1]
    return trend



def find_valid_domain(da):
    """ finds first and last index of non-NaN values for zonal integrals 
    
    input:
    da                       .. xr DataArray
    
    output:
    min/max_lat (/max_depth) .. indices of first/last non-NaN entry in da 
    """
        
    if da.ndim==2:  # lat_bin field
        min_lat = 0 
        while np.isnan(da[0,min_lat]).item():
            min_lat += 1

        max_lat = len(da[0,:])-1
        while np.isnan(da[0,max_lat]).item():
            max_lat -= 1
        print(min_lat, max_lat)
        return (min_lat, max_lat)
    
    if da.ndim==3:  # depth-lat_bin field
        min_lat = 0 
        while np.isnan(da[0,0,min_lat]).item():
            min_lat += 1

        max_lat = len(da[0,0,:])-1
        while np.isnan(da[0,0,max_lat]).item():
            max_lat -= 1
        
        
        max_depth = 41
        while np.isnan(da[0,max_depth, min_lat]).item():
            max_depth -= 1
        return (min_lat, max_lat, max_depth)