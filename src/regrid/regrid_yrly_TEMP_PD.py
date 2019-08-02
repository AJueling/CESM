import os
import sys
sys.path.append("..")
import numpy as np
import xarray as xr

from paths import path_samoc

def add_attrs(ds):
    ds.TEMP.attrs = {'long_name':     'Potential Temperature',
                     'units':         'degC',
                     'grid_loc':      '3111',
                     'cell_methods':  'time: mean'}
    ds.PD.attrs   = {'long_name':     'Potential Density Ref to Surface',
                     'units':         'gram/centimeter^3',
                     'grid_loc':      '3111',
                     'cell_methods':  'time: mean'}
    return ds


for y in np.arange(71,301):
    print(y)
    fn_orig = f'{path_samoc}/ctrl/ocn_yrly_TEMP_PD_{y:04d}.nc'
    fn_temp = f'{path_samoc}/ctrl_rect/TEMP_PD_yrly_{y:04d}.nc'
    assert os.path.exists(fn_orig)
    ds = xr.open_dataset(fn_orig, decode_times=False)
    ds = add_attrs(ds)
    ds.to_netcdf(fn_temp)
    cmd = f'/home/dijkbio/andre/tx0.1_to_0.4/nccurv2ncrect.sc {fn_temp} {path_samoc}/ctrl_rect'
    
    os.system(cmd)
    os.remove(fn_temp)