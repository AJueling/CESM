""" creating CICE_XMXL plots """

import sys
sys.path.append("..")
run     = sys.argv[1]

import xarray as xr
import datetime

from grid import generate_lats_lons
from CICE import CICE_XMXL_plots
from paths import CESM_filename, file_ex_ocn_rcp
from timeseries import IterateOutputCESM

if __name__=="__main__":
    print(f'running CICE_XMXL plots run={run}')
    print(f'{datetime.datetime.now()}\n\n')

    lats, lons = generate_lats_lons('ocn')
    MASK = xr.open_dataset(file_ex_ocn_rcp, decode_times=False).REGION_MASK
    
    for i, (y,m,ocn_file) in enumerate(IterateOutputCESM(domain='ocn', run=run, tavg='monthly')):
        if m==1: print(y)
        ice_file = CESM_filename(domain='ice', run=run, y=y, m=m)
        aice = xr.open_dataset(ice_file, decode_times=False).aice[0,:,:]
        XMXL = xr.open_dataset(ocn_file, decode_times=False).XMXL[0,:,:]
        CICE_XMXL_plots(aice=aice, XMXL=XMXL, lons=lons, lats=lats, MASK=MASK, run=run, i=i)
        aice.close()
        XMXL.close()
        
        
    print(f'\n\nfinished at\n{datetime.datetime.now()}')