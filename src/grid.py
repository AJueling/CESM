import numpy as np
import xarray as xr
from operations.read_binary import read_binary_2D_double

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
    

def generate_lats_lons_CESM(grid_file):
    """
    genrates lats and lons fields (no shift)
    """
    imt,jmt = 3600,2400
    lats = read_binary_2D_double(grid_file,imt,jmt,1)
    lons = read_binary_2D_double(grid_file,imt,jmt,2)
    
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