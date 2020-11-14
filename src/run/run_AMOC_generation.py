import sys
import datetime
sys.path.append("..")
import numpy as np
import xarray as xr

from MOC import calculate_AMOC_sigma_z
import tqdm.notebook as tqdm_notebook

from paths import path_results, path_prace, file_RMASK_ocn, file_RMASK_ocn_low, file_ex_ocn_ctrl, file_ex_ocn_lpd

if __name__=='__main__':
    """ calculate the AMOC streamfunctions """
    print(f'{datetime.datetime.now()}\n\n')
    for run in ['lpd', 'lr1', 'ctrl', 'rcp']:
        if run in ['ctrl', 'rcp']:   domain, fe = 'ocn', file_ex_ocn_ctrl
        elif run in ['lpd', 'lr1']:  domain, fe = 'ocn_low', file_ex_ocn_lpd
        if run=='ctrl':  yy = [227,228,229]# np.arange(200, 230)
#         elif run=='lpd':  yy = np.arange(500, 530)
        if run in ['rcp']:   yy = [2029,2030,2031,2032,2033,2093,2094,2095,2096,2097,2098,2099,2100,2101,2062,2063,2064,2065,2066]# np.arange(2067, 2101)  # 
        else: continue
        geometry = xr.open_dataset(fe, decode_times=False)[['DXT', 'DYT', 'DXU', 'DYU', 'REGION_MASK']].squeeze()
        for y in tqdm_notebook.tqdm(yy):
            print(run, y)
            fn = f'{path_prace}/MOC/AMOC_sz_yz_{run}_{y}.nc'
    #         if os.path.exists(fn):  continue
            PD = xr.open_dataset(f'{path_prace}/{run}/ocn_yrly_TEMP_PD_{y:04d}.nc', decode_times=False)['PD'].drop(['TLAT', 'TLONG', 'ULAT', 'ULONG'])
            VVEL = xr.open_dataset(f'{path_prace}/{run}/ocn_yrly_UVEL_VVEL_{y:04d}.nc', decode_times=False)['VVEL'].drop(['TLAT', 'TLONG', 'ULAT', 'ULONG'])
            ds = xr.merge([geometry, PD, VVEL])
            AMOC_yz, AMOC_sz = calculate_AMOC_sigma_z(domain=domain, ds=ds, fn=fn)
            
    print(f'\n\nfinished at\n{datetime.datetime.now()}')    