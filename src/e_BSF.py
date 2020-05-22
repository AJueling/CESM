import sys
import numpy as np
import xarray as xr

from paths import file_ex_ocn_ctrl, file_ex_ocn_lpd
from regions import Drake_Passage, DP_North, WG_center
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


if __name__=='__main__':

    run = sys.argv[1]
    # ys, ye = int(sys.argv[2]), int(sys.argv[3])
    if run=='ctrl':              yy=np.arange(200,230)
    elif run=='lpd':             yy=np.arange(500,530)
    elif run in ['rcp', 'lr1']:  yy = np.arange(2000,2101)

    if run in ['ctrl', 'rcp']:   
        DZU = xr_DZ('ocn', grid='U')
        DYU = xr.open_dataset(file_ex_ocn_ctrl, decode_times=False).DYU
    elif run in ['lpd', 'lr1']:  
        DZU = xr_DZ('ocn_low', grid='U')
        DYU = xr.open_dataset(file_ex_ocn_lpd, decode_times=False).DYU
    for i, y in enumerate(yy):
        UVEL = xr.open_dataset(f'{path_prace}/{run}/ocn_yrly_UVEL_VVEL_{y:04d}.nc', decode_times=False).UVEL

        psi_list.append( (DYU*(UVEL*DZU).sum('z_t')/1e10).cumsum('nlat'))
    psi = xr.concat(psi_list, concat_dim='time')
    psi.to_netcdf(f'{path_prace}/BSF/BSF_{run}.nc')