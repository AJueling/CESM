import numpy as np
import xarray as xr

from regions import Drake_Passage


def calculate_BSF(ds):
    """ barotropic stream function [m^3/s] 
    
    input:
    ds  .. xr Dataset
    """
    assert 'DYU' in ds
    assert 'UVEL' in ds
    assert 'REGION_MASK' in ds
    
    DZU = xr_DZ('ocn', grid='U')
    BSF = (ds.DYU*ds.UVEL[0,:,:,:]*DZU).where(ds.REGION_MASK>0).sum(dim='z_t')/1e4
    for j in np.arange(1,2400):
        BSF[j,:] += BSF[j-1,:]
        
    return BSF


def DP_transport(BSF):
    """ Drake Passage transport  [m^3/s]
    
    input:
    BSF .. 
    """
    DPT = BSF.sel(Drake_Passage)[-1] - BSF.sel(Drake_Passage)[0]
    return DPT