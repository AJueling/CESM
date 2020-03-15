import os
import sys
import scipy as sp
import numpy as np
import pandas as pd
import xarray as xr
import cmocean
import cartopy
import cartopy.crs as ccrs
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('Agg')

start = int(sys.argv[1])
end   = int(sys.argv[2])

sys.path.append("..")
from paths import path_data, path_results
from paths import file_ex_ocn_ctrl, file_ex_ocn_lpd
from paths import file_RMASK_ocn, file_RMASK_ocn_rect
from regions import AMO_mask, boolean_mask

ctrl_daily = xr.open_mfdataset('/projects/0/prace_imau/prace_2013081679/cesm1_0_4/spinup_pd_maxcores_f05_t12/OUTPUT/ocn/hist/daily/spinup_pd_maxcores_f05_t12.pop.h.nday1.0300*.nc', combine='nested', concat_dim='time').SST

MASK = boolean_mask(domain='ocn', mask_nr=0) + boolean_mask(domain='ocn', mask_nr=-13) + boolean_mask(domain='ocn', mask_nr=-14)

cmap = plt.cm.Spectral_r
cmap.set_under('w')
daterange = pd.date_range('2018-01-01', '2018-12-31')
for i in range(365):
#     fn =f'{path_results}/NAC_presentation/global_daily_SST_ctrl_{i:03d}.png'
    fn =f'{path_results}/NAC_presentation/NA_daily_SST/NA_daily_SST_ctrl_black_{i:03d}.png'
    if i not in np.arange(start,end):  continue
    if os.path.exists(fn):  continue
#     fig = plt.figure(figsize=(15, 7))
#     ax = plt.axes(projection=ccrs.PlateCarree())
    
    fig = plt.figure(figsize=(7, 7))
    fig.patch.set_facecolor('k')
    ax = plt.axes(projection=ccrs.NearsidePerspective(central_longitude=-40.0, central_latitude=40.0))
    
    ax.set_position((0,0,1,1))
    im = ax.pcolormesh(ctrl_daily.TLONG, ctrl_daily.TLAT, ctrl_daily.isel(time=i).where(MASK),
                       cmap=cmap, vmin=-1.7, vmax=31, transform=ccrs.PlateCarree(),
                       )
    ax.add_feature(cartopy.feature.LAND, zorder=2, edgecolor='black', facecolor='grey')

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
    gl.n_steps = 90
    gl.ylocator = matplotlib.ticker.FixedLocator(np.arange(-90,100,30))
    gl.xlocator = matplotlib.ticker.FixedLocator(np.arange(0,370,30))
    fig.text(.02, .94, daterange[i].strftime('%b %d'), fontsize=14, color='grey')
    plt.savefig(fn, dpi=150, transparent=True)
    plt.close()