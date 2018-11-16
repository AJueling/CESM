import sys
import datetime
sys.path.append("..")
import xarray as xr

from OHC import OHC_vert_diff_mean_rm
from paths import path_samoc

run     = sys.argv[1]
domain  = 'ocn'

assert run in ['ctrl', 'rcp']

print(f'running OHC_vert_diff_mean_rm.py run={run}')
print(f'{datetime.datetime.now()}\n\n')

if run=='ctrl':
    ds = xr.open_dataset(f'{path_samoc}/OHC/OHC_integrals_ctrl_Global_Ocean.nc'  , decode_times=False)
    for field in ['OHC_global', 'OHC_global_levels', 'OHC_zonal', 'OHC_zonal_levels']:
        ds[field][105] = (ds[field].sel({'time':204*365}) + ds[field].sel({'time':206*365}) )/2
elif run=='rcp':
    ds = xr.open_dataset(f'{path_samoc}/OHC/OHC_integrals_rcp_Global_Ocean.nc'  , decode_times=False)
    
OHC_vert_diff_mean_rm(ds=ds, run=run)

print(f'\n\nfinished at\n{datetime.datetime.now()}')