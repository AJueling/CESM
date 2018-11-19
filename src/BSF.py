import numpy as np
import xarray as xr

from regions import Drake_Passage, DP_North, WG_center
# from timeseries import round_tlat_tlong
from xr_DataArrays import xr_DYU, xr_DZ


def calculate_BSF(ds):
    """ barotropic stream function [m^3/s] 
    
    input:
    ds  .. xr Dataset
    """
    assert 'UVEL' in ds
    
    DYU = xr_DYU('ocn')
    # remocing unnecessary TLAT/TLONG because they are not always exactly equal
    # (due to rounding errors) and hence lead to Error messages
    DYU = DYU.drop(['TLAT', 'TLONG'])
    
    DZU = xr_DZ('ocn', grid='U')
    
    UVEL = ds.UVEL/1e2  # [cm] to [m]
    if 'TLAT'  in UVEL.coords:  UVEL = UVEL.drop('TLAT')
    if 'TLONG' in UVEL.coords:  UVEL = UVEL.drop('TLONG')
    
    BSF = (((UVEL*DZU).sum(dim='z_t'))*DYU)
    for j in np.arange(1,2400):
        BSF[j,:] += BSF[j-1,:]
        
    return BSF


def DP_transport(BSF):
    """ Drake Passage transport  [m^3/s]
    
    input:
    BSF .. xr DataArray
    """
    DPT = BSF.sel(DP_North)
    return DPT


def WG_transport(BSF):
    """ Weddell Gyre transport [m^3/s]
    
    input:
    BSF .. xr DataArray
    """
    WGT = BSF.sel(WG_center)
    return WGT