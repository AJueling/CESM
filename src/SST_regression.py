""" SST regression pattern analysis

this script analyzes the regression patterns obtained by `SST_generation.py` and

1. plots of regression pattern
2. quantitative stationarity check of the regression patterns
3. comparison of the simulated to the observed pattern (with spatial correlations)

"""

import os
import sys
import numpy as np
import xarray as xr
import subprocess
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('agg')
# matplotlib.rc_file('rc_file')

from maps import regr_map
from paths import path_samoc, path_results, file_ex_ocn_ctrl
from bc_analysis_fields import AnalyzeField as AF

run = str(sys.argv[1])

assert run in ['ctrl', 'lpd']

if run=='ctrl':   starts = np.arange(50, 151, 5)
elif run=='lpd':  starts = np.arange(154, 415, 10)

# ==================================================================================================
print(f'Analysis of regression patterns for {run}')
# ==================================================================================================

if run=='ctrl':
    TLAT = xr.open_dataset(file_ex_ocn_ctrl, decode_times=False).TLAT.coords['TLAT']

for i, idx in enumerate(['AMO', 'SOM', 'TPI']):
    had   = xr.open_dataset(f'{path_samoc}/SST/{idx}_regr_had.nc').slope
    region = [{'longitude':slice(-80,0), 'latitude':slice(60,0)}, 1, 2][i]
    
    fn_new = f'{path_samoc}/SST/{idx}_spatial_correlations_{run}.nc'
    da = xr.DataArray(data=np.zeros(len(starts)),
                      coords={'time': starts},
                      dims=('time'))
    
    for k, t in enumerate(starts):
        tslice = (t, t+148)
    
        try:
            fn = f'{path_samoc}/SST/{idx}_regr_{run}_{tslice[0]}_{tslice[1]}.nc'
            assert os.path.exists(fn)
            ds = xr.open_dataset(fn, decode_times=False)
            
        except:
            print(f'{idx} file for segment {tslice[0]}-{tslice[1]} does not exist.')
            continue
        
        # visualization
        print('plotting')
        fn_map = f'{path_results}/SST/{idx}_regr_map_{run}_{tslice[0]}_{tslice[1]}'
        try:
            assert os.path.exists(fn_map)
        except:
            regr_map(ds=ds, index=idx, run=run, fn=fn_map)
    
        # compare to observations
        print('calculating spatial correlations')
        
        segment = ds.slope
        if run=='ctrl':
            segment.coords['TLAT'] = TLAT
            
        da.values[k] = AF().spatial_correlation(field_A=had, field_B=segment,
                                                    selection=region)
    da.to_netcdf(fn_new)
    
    # spatial correlation plot
    plt.figure()
    da.plot()
    plt.savefig(f'{path_results}/SST/{idx}_spatial_correlation(t)_{run}')
        
    # video of maps
    fn_temp = 'run.sh'
    with open(fn_temp, 'w') as file:
        file.write('#!/bin/bash\n')
        file.write('module load ffmpeg\n')
        file.write('cd /home/dijkbio/andre/CESM/results/SST\n')
        file.write(f'ffmpeg -i {idx}_regr_map_{run}_%*.png -framerate 2 {idx}_regr_maps_{run}.gif\n')
        file.write(f'ffmpeg -i {idx}_regr_map_{run}_%*.png -framerate 2 -qscale:v 0 {idx}_regr_maps_{run}.avi')
        file.close()
    os.chmod(fn_temp, 0o755)
    subprocess.call("./"+fn_temp)
    os.remove(fn_temp)