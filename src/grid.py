import os
import numpy as np
import xarray as xr

from paths import path_results, file_ex_ocn_ctrl, grid_file
from constants import imt, jmt
from read_binary import read_binary_2D_double
from xr_DataArrays import xr_DZ

def generate_lats_lons(grid_file):
    """
    genrates lats and lons fields and shifts them so they are increasing 
    (important for plotting with Basemap)
    """
    imt,jmt = 3600,2400
    lats = read_binary_2D_double(grid_file,imt,jmt,1)
    lons = read_binary_2D_double(grid_file,imt,jmt,2)
    
    shift = np.zeros((jmt),dtype='int')
    for j in range(jmt):
        if j<jmt-1:  b = imt-np.argmin(lons[:,j])
        if j==jmt-1: b = 900
        lats[:,j] = 180/np.pi*np.roll(lats[:,j],b)
        lons[:,j] = 180/np.pi*np.roll(lons[:,j],b)
        shift[j]  = b
    lons[imt-1,jmt-1] = 90.
    
    return lats, lons, shift
    

def generate_lats_lons(domain):
    """
    genrates lats and lons fields (no shift)
    """
    assert domain in ['ocn']
    lats = read_binary_2D_double(grid_file,imt,jmt,1)
    lons = read_binary_2D_double(grid_file,imt,jmt,2)
    lats = lats.T*180/np.pi  # rad to deg
    lons = np.roll(lons*180/np.pi+180, int(imt/2), axis=0).T  # rad to deg; same 
    return lats, lons 


def shift_field(field,shift):
    """
    shifts a 2D (imt,jmt) field
    """
    imt,jmt = 3600,2400
    shifted = np.zeros((imt,jmt))
    for j in range(jmt):
        shifted[:,j]  = np.roll(field[:,j],shift[j])
    return shifted


def create_dz_mean(domain):
    """ average depth [m] per level of """
    assert domain=='ocn'
    
    fn = f'{path_results}/geometry/dz_mean'
    if os.path.exists(fn):
        dz_mean = xr.open_dataarray(fn, decode_times=False)
    else:
        DZT     = xr_DZ(domain)
        dz_mean = DZT.where(DZT>0).mean(dim=('nlat', 'nlon'))
        dz_mean.to_netcdf(fn)
        
    return dz_mean


def create_tdepth(domain):
    """
    input:
    domain .. 'ocn'
    
    output:
    tdepth .. np array
    """
    assert domain=='ocn'
    
    fn = f'{path_results}/geometry/tdepth.csv'
    if os.path.exists(fn):
        tdepth = np.genfromtxt(fn, delimiter=',')
    else:
        ds = xr.open_dataset(file_ex_ocn_ctrl, decode_times=False)
        tdepth = ds.z_t.values/1e2
        np.savetxt(fn, tdepth, delimiter=',')
        
    return tdepth


def find_array_idx(array, val):
    """ index of nearest value in array to val 
    
    input:
    array .. array like
    value .. value to be approximately in array
    """
    idx = (np.abs(array - val)).argmin()
    return idx
    