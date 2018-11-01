import numpy as np
import xarray as xr

from regions import Drake_Passage, boolean_mask, DP_transport
from xr_DataArrays import xr_DYU, xr_DZ


def calculate_BSF(ds):
    """ barotropic stream function [m^3/s] 
    
    input:
    ds  .. xr Dataset
    """
    assert 'UVEL' in ds
    
    DYU = xr_DYU('ocn')
    DZU = xr_DZ('ocn', grid='U')
    GLOBAL_MASK = boolean_mask('ocn', 0)
    
    BSF = (((ds.UVEL*DZU).sum(dim='z_t'))*DYU).where(GLOBAL_MASK)/1e4
    for j in np.arange(1,2400):
        BSF[j,:] += BSF[j-1,:]
        
    return BSF


def DP_transport(BSF):
    """ Drake Passage transport  [m^3/s]
    
    input:
    BSF .. 
    """
    DPT = BSF.sel(DP_transport) #- BSF.sel(Drake_Passage)[0]
    return DPT