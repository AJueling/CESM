import os
import sys
import dask
import numpy as np
import scipy as sp
import xarray as xr

sys.path.append("..")
from tqdm import tqdm_notebook
from LFCA import xlfca
from paths import path_prace, path_data, path_results
from filters import lowpass, Lanczos, deseasonalize
from regions import boolean_mask, mask_box_in_region
from eofs.xarray import Eof
from xr_DataArrays import xr_AREA
from ab_derivation_SST import DeriveSST as DS
from bb_analysis_timeseries import AnalyzeTimeSeries as ATS

run   = sys.argv[1]
basin = sys.argv[2]
dsdt  = sys.argv[3]

if run=='had':     dt, domain = ''        , 'ocn_had'
elif run=='ctrl':  dt, domain = '_51_301' , 'ocn_rect'
elif run=='lpd':   dt, domain = '_154_404', 'ocn_low'
    
if basin=='North_Pacific':   mask_nr, bounding_lats, bounding_lons = 2, ( 20,68), (110,255)
if basin=='full_Pacific':    mask_nr, bounding_lats, bounding_lons = 2, (-38,68), (110,290)
if basin=='Southern_Ocean':  mask_nr, bounding_lats, bounding_lons = 1, None    , None
if basin=='North_Atlantic':  mask_nr, bounding_lats, bounding_lons = 6, (  0,60), (-80,  0)
fn = f'{path_prace}/SST/SST_monthly_{dsdt}_{run}{dt}.nc'

MASK = mask_box_in_region(domain=domain, mask_nr=mask_nr, bounding_lats=bounding_lats, bounding_lons=bounding_lons)
AREA = xr_AREA(domain=domain).where(MASK)
SST = xr.open_dataarray(fn, decode_times=False).where(MASK)
if basin in ['North_Pacific', 'full_Pacific'] and run=='had':  # shifting
    AREA = DS().shift_had(AREA)
    SST = DS().shift_had(SST)
if basin=='North_Atlantic' and run=='ctrl':
    AREA = DS().shift_ocn_rect(AREA)
    SST = DS().shift_ocn_rect(SST)
if basin=='North_Atlantic' and run=='lpd':
    AREA = DS().shift_ocn_low(AREA)
    SST = DS().shift_ocn_low(SST)
AREA = AREA.where(np.isnan(SST[0,:,:])==False, drop=True)
SST = SST.where(np.isnan(SST[0,:,:])==False, drop=True)
#             AREA.plot(ax=ax[k,i])
scale = AREA/AREA.sum()
scale = xr.apply_ufunc(np.sqrt, scale)
for n_EOFs in [3, 30]:
    fn_lfca = f'{path_prace}/LFCA/LFCA_{run}_{basin}_{dsdt}_n{n_EOFs}.nc'
    if os.path.exists(fn_lfca):  continue
    lfca = xlfca(x=SST, cutoff=120, truncation=n_EOFs, scale=scale)
    lfca.to_netcdf(fn_lfca)