import sys
import numpy as np
import xarray as xr
import datetime

from eofs.xarray import Eof

from paths import path_samoc, path_results

from timeseries import IterateOutputCESM
from regions import boolean_mask, TPI_masks
from xr_DataArrays import xr_AREA

# sys.path.append("..")



def SST_index(xa_SST, AREA, index_loc, AREA_index, MASK, dims=('nlat', 'nlon')):
    """ calculates the average SST over an area, possibly as a time series """
    assert type(xa_SST)==xr.core.dataarray.DataArray
    assert type(AREA)==xr.core.dataarray.DataArray
    if type(index_loc)==dict:
        index = (xa_SST*AREA).where(MASK).sel(index_loc).sum(dim=dims)/AREA_index
    if index_loc==None:
        index = (xa_SST*AREA).where(MASK).sum(dim=dims)/AREA_index
    return index


def SST_index_from_monthly(run, index_loc, MASK):
    """ loads monthly SST data, calculated SST_index, returns raw timeseries"""
    assert run in ['ctrl', 'rcp', 'lpd', 'lpi']
    if run in ['ctrl', 'rcp']:
        domain = 'ocn'
    elif run in ['lpd', 'lpi']:
        domain = 'ocn_low'
    AREA = xr_AREA(domain)
    AREA_index = AREA.where(MASK).sum()
    
    for  i, (y, m, s) in enumerate(IterateOutputCESM(domain=domain, run=run, tavg='monthly')):
        if m==1: print(y)
        xa_SST = xr.open_dataset(s, decode_times=False).TEMP[0,0,:,:]
        SSTi = SST_index(xa_SST, AREA, index_loc, AREA_index, MASK)
        if i==0:
            new_SSTi = SSTi.copy()
        else:
            new_SSTi = xr.concat([new_SSTi,SSTi], dim='time')
            
    return new_SSTi


def EOF_SST_analysis(xa, weights, neofs=1, npcs=1, fn=None):
    assert type(xa)==xr.core.dataarray.DataArray
    assert type(weights)==xr.core.dataarray.DataArray
    assert 'time' in xa.dims
    assert np.shape(xa[0,:,:])==np.shape(weights)
    
    # anomalies by removing time mean
    xa = xa - xa.mean(dim='time')
    # Retrieve the leading EOF, expressed as the covariance between the leading PC
    # time series and the input xa anomalies at each grid point.
    solver = Eof(xa, weights=weights)
    eofs = solver.eofsAsCovariance(neofs=neofs)
    pcs = solver.pcs(npcs=npcs, pcscaling=1)
    if fn!=None:
        xr.merge([eofs, pcs]).to_netcdf(fn)
    return eofs, pcs


def IPO_TPI(run):
    """ raw timeseries of the three tripole SST indices described in Henley et al. (2015) """
    assert run in ['ctrl', 'rcp', 'lpd', 'lpi']
    if run in ['ctrl', 'rcp']:
        domain = 'ocn'
    elif run in ['lpd', 'lpi']:
        domain = 'ocn_low'
    
    TPI1 = SST_index_from_monthly(run=run, index_loc=None, MASK=TPI_masks(domain, 1))
    TPI2 = SST_index_from_monthly(run=run, index_loc=None, MASK=TPI_masks(domain, 2))
    TPI3 = SST_index_from_monthly(run=run, index_loc=None, MASK=TPI_masks(domain, 3))
    
    TPI1.to_netcdf(f'{path_samoc}/SST/TPI1_{run}.nc')
    TPI2.to_netcdf(f'{path_samoc}/SST/TPI2_{run}.nc')
    TPI3.to_netcdf(f'{path_samoc}/SST/TPI3_{run}.nc')
    
    return
    

if __name__=="__main__":
    idx = sys.argv[1]
    run = sys.argv[2]
    
    assert idx in ['IPO']
    assert run in ['ctrl', 'rcp', 'lpd', 'lpi']
    
    print(f'running {idx} index of {run} run\n')
    print(f'started at\n{datetime.datetime.now()}\n\n')

    if idx=='IPO':
        IPO_TPI(run)
        
    print(f'\n\nfinished at\n{datetime.datetime.now()}')