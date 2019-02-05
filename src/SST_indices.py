
import xarray as xr

from SST import SST_index
from paths import path_samoc, path_results
from regions import mask_box_in_region, bll_SOM
from timeseries import deseasonalize
from xr_DataArrays import xr_AREA

"""
goal: perform EOF of monthly, deseasonalized SST (detrended with global mean SST)
NB: we will use 60S to 60N SST data as polar data is limited in observations

1. for detrending: compute monthly global mean SST time series, deseasonalize them
    >>> see `SST_data_generation.py` file for scripts
2. North Pacific monthly output fields
    2.1. create monthly SST field
         (if appropriate: determine extend of grid, limit all coordinates)
         save as single file
         a) North of 38 deg S
         b) North of Equator
         b) North of 20 deg N
    2.2. deseasonalize 
    2.3. detrend global mean, deseasonalized SST, 
    2.4. (remove mean at each point)
3. EOF analysis
   --> index is first principal component
4. regress time series on global maps
"""

def PMV_EOF_indices(run, extent):
    """
    
    """
    assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'had']
    assert extent in ['38S', 'Eq', '20N']
    if run in ['ctrl', 'rcp']:
        domain = 'ocn_rect'
        run_name = f'rect_{run}'
    elif run in ['lpd', 'lpi']:
        domain = 'ocn_low'
        run_name = run
    elif run=='had':
        domain = 'ocn_had'
        run_name = run
    
    # load fields (generated in `SST_data_generatoin.py`)
    # <SST>(t):  60S-60N (tropical/extratropical) monthly timeseries
    SST_xm = xr.open_dataarray(f'{path_samoc}/SST/SST_60S_60N_mean_monthly_{run_name}.nc', decode_times=False)
    
    # SST(t,x,y): monthly SST field limited to Pacific North of `extent` (*)
    SST = xr.open_dataarray(f'{path_samoc}/SST/SST_monthly_{run_name}.nc', decode_times=False)
    
    # deseasonalize
    SST_xm_ds = deseasonalize(SST_xm )
    SST_ds    = deseasonalize(SST)
    
    # some of the time series are not the same length
    if run=='ctrl' :  SST_xm_ds = SST_xm_ds[:-7]
    elif run=='rcp':  SST_xm_ds = SST_xm_ds[:-1]
    
    # detrend
    SST_ds_dt = SST_ds - SST_xm_ds
    
    # remove mean at each point
    SST_ds_dt_dm = SST_ds_dt - SST_ds_dt.mean('time')
    
    # EOF analysis
    # N.B. cut off two year on either end as the arbitrary start month biases the filtered time series
    if extent=='38S':
        latS, lonE = -38, 300
    elif extent=='Eq':
        latS, lonE = 0, 285
    elif extent=='20N':
        latS, lonE = 20, 255
    
    AREA = xr_AREA(domain=domain)
    Pac_MASK = mask_box_in_region(domain=domain, mask_nr=2, bounding_lats=(latS,68), bounding_lons=(110,lonE))
    Pac_area = AREA.where(Pac_MASK)
    fn = f'{path_samoc}/SST/SST_PDO_EOF_{extent}_{run}.nc'
    
    eof, pc  = EOF_SST_analysis(xa=SST_ds_dt_dm[24:-24,:,:].where(Pac_MASK),
                                weights=Pac_area, fn=fn)
    
    return eof, pc



def SOM_index(run):
    """ """
    assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'had']
    print(run)
    
    if run in ['ctrl', 'rcp']:
        domain = 'ocn'
        dims = ('nlat', 'nlon')
    elif run in ['lpd', 'lpi']:
        domain = 'ocn_low'
        dims = ('nlat', 'nlon')
    elif run=='had':
        domain = 'ocn_had'
        dims = ('latitude', 'longitude')
        
    (blats, blons) = bll_SOM
    
    # load yearly data files
    
    MASK = mask_box_in_region(domain=domain, mask_nr=0, bounding_lats=blats, bounding_lons=blons)
    AREA = xr_AREA(domain=domain).where(MASK)
    SOM_area = AREA.sum()
    
    SST_yrly = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_{run}.nc').where(MASK)
    SOM = SST_index(xa_SST=SST_yrly, AREA=AREA, AREA_index=SOM_area, MASK=MASK, dims=dims)
    SOM.to_netcdf(f'{path_samoc}/SST/SOM_raw_{run}.nc')
    
    return SOM
