import sys
import datetime
sys.path.append("..")
import xarray as xr

from OHC import OHC_vert_diff_mean_rm

run     = sys.argv[1]
domain  = 'ocn'

assert run in ['ctrl', 'rcp']

print('run_OHC_vert_diff_mean_rm.py')
print(datetime.datetime.now())

if run=='ctrl':
    ds = xr.open_dataset(f'{path_samoc}/OHC/OHC_integrals_ctrl_Global_Ocean.nc'  , decode_times=False)
    for field in ['OHC_global', 'OHC_global_levels', 'OHC_zonal', 'OHC_zonal_levels']:
        ds[field][105] = (ds[field].sel({'time':204*365}) + ds[field].sel({'time':206*365}) )/2
else:
    ds = xr.open_dataset(f'{path_samoc}/OHC/OHC_integrals_rcp_Global_Ocean.nc'  , decode_times=False)
    
OHC_vert_diff_mean_rm(ds=ds, run=run)

print(datetime.datetime.now())