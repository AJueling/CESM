import numpy as np
import xarray as xr

from grid import create_tdepth, find_array_idx
from paths import file_ex_ocn_ctrl
from xr_DataArrays import xr_DZ, xr_DXU


def calculate_MOC(ds, DXU, DZU, MASK):
    """ Atlantic Meridional Overturning circulation 
    
    input:
    ds   .. xr Dataset of CESM output
    
    output:
    MOC .. 2D xr DataArray
    """
    assert 'VVEL' in ds
    MOC = (ds.VVEL*DXU*DZU).where(MASK).sum(dim='nlon')/1e2  # [m^3/s]
    for k in np.arange(1,42):
        MOC[k,:] += MOC[k-1,:]
    
    return MOC


def approx_lats(domain):
    """ array of approx. latitudes for ocean """
    assert domain=='ocn'
    ds = xr.open_dataset(file_ex_ocn_ctrl, decode_times=False)
    lats = ds.TLAT[:,900].copy()
    lats[:120] = -78
    return lats



def AMOC_max(AMOC):
    """ AMOC maximum at 26 deg N, 1000 m """
    lats    = approx_lats('ocn')
    tdepths = create_tdepth('ocn')
    j26   = find_array_idx(lats, 26)
    z1000 = find_array_idx(tdepths, 1000)
    
    return AMOC.isel({'z_t':z1000, 'nlat':j26})



def plot_AMOC(AMOC):
    lats    = approx_lats('ocn')
    tdepths = create_tdepth('ocn')
    j26   = find_array_idx(lats, 26)
    z1000 = find_array_idx(tdepths, 1000)
    return