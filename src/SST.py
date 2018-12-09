import numpy as np
import xarray as xr
from eofs.xarray import Eof

def SST_index(xa_SST, AREA, index_loc, AREA_index, MASK, dims=('nlat', 'nlon')):
    """ calculates the average SST over an area, possibly as a time series """
    assert type(xa_SST)==xr.core.dataarray.DataArray
    assert type(AREA)==xr.core.dataarray.DataArray
    if type(index_loc)==dict:
        index = (xa_SST*AREA).where(MASK).sel(index_loc).sum(dim=dims)/AREA_index
    if index_loc==None:
        index = (xa_SST*AREA).where(MASK).sum(dim=dims)/AREA_index
    return index


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