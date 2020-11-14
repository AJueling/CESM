import sys
import datetime
sys.path.append("..")
import xarray as xr
from paths import path_prace
from regions import boolean_mask, Atlantic_mask
from timeseries import IterateOutputCESM
from xr_DataArrays import xr_AREA


for j, run in enumerate(['ctrl','lc1']):
    # if j==1: break
    domain = ['ocn', 'ocn_low'][j]
    TAREA = xr_AREA(domain=domain)
    mask_A = Atlantic_mask(domain=domain)
    mask_P = boolean_mask(domain=domain, mask_nr=2)
    mask_S = boolean_mask(domain=domain, mask_nr=1)
    shf, shf_A, shf_P, shf_S = [], [], [], []
    for i, (y,m,f) in enumerate(IterateOutputCESM(domain=domain, run=run, tavg='monthly')):
        da = xr.open_dataset(f, decode_times=False).SHF*TAREA
        shf.append(da.sum())
        shf_A.append(da.where(mask_A).sum())
        shf_P.append(da.where(mask_P).sum())
        shf_S.append(da.where(mask_S).sum())
        # if i==24:  break
    shf = xr.concat(shf, dim='time')
    shf_A = xr.concat(shf_A, dim='time')
    shf_P = xr.concat(shf_P, dim='time')
    shf_S = xr.concat(shf_S, dim='time')
    shf.name = 'Global_Ocean'
    shf_A.name = 'Atlantic_Ocean'
    shf_P.name = 'Pacific_Ocean'
    shf_S.name = 'Southern_Ocean'
    xr.merge([shf, shf_A, shf_P, shf_S]).to_netcdf(f'{path_prace}/SHF/SHF_monthly_{run}.nc')
