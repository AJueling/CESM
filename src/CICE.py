import os
import sys
sys.path.append("..")
from shutil import copyfile

import numpy as np
import xarray as xr
import cmocean
import cartopy
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from paths import path_results



def CICE_XMXL_plots(aice, XMXL, lons, lats, MASK, run, i):
    """"""
    year  = int(XMXL.time.item()/365)
    month = int((XMXL.time.item()%365)/30)
    filename = f'{path_results}/CICE/CICE_XMXL_video/CICE_XMXL_{run}_{year}_{month}.png'
    if os.path.exists(filename):
        print(f'{filename} exists already')
        return

    f   = plt.figure(figsize=(8,8))
    ax1 = plt.subplot(2,2,1, projection=ccrs.NearsidePerspective(                      central_latitude= 90))
    ax2 = plt.subplot(2,2,2, projection=ccrs.NearsidePerspective(                      central_latitude=-90))
    ax3 = plt.subplot(2,2,3, projection=ccrs.NearsidePerspective(satellite_height=1e7, central_latitude= 90))
    ax4 = plt.subplot(2,2,4, projection=ccrs.NearsidePerspective(satellite_height=1e7, central_latitude=-90))


    for ax in [ax1, ax2]:  # MXL
        im = ax.pcolormesh(XMXL.TLONG, XMXL.TLAT, XMXL.where(MASK>0)/100,
                           cmap=cmocean.cm.amp, vmin=0, vmax=500,
                           transform=ccrs.PlateCarree() )

    cbar_ax1 = f.add_axes([0.48, 0.57, 0.02, 0.35])
    f.colorbar(im, cax=cbar_ax1, extend='max')

    for ax in [ax3, ax4]:  # ICE
        im = ax.pcolormesh(lons, lats, aice,
                           cmap=cmocean.cm.ice, vmin=0, vmax=100,
                           transform=ccrs.PlateCarree() )

    cbar_ax2 = f.add_axes([0.48, 0.08, 0.02, 0.35])
    f.colorbar(im, cax=cbar_ax2)

    cbar_ax1.text(.5, 1.08, f'{year} - {month:2d}'         , ha='center', transform=cbar_ax1.transAxes, fontsize=18)
    cbar_ax1.text(.5, -.16, 'maximum mixed layer depth [m]', ha='center', transform=cbar_ax1.transAxes, fontsize=14)
    cbar_ax2.text(.5, 1.06, 'sea ice cover [%]'            , ha='center', transform=cbar_ax2.transAxes, fontsize=14)

    for ax in [ax1, ax2, ax3, ax4]:
        ax.add_feature(cartopy.feature.LAND, zorder=2, edgecolor='black', facecolor='lightgrey')
    #     ax.coastlines(resolution='110m')
        ax.gridlines()

    plt.subplots_adjust(left=0.02, right=0.98, top=0.95, bottom=0.05)
    
    f.savefig(filename)
    
    copyfile(filename, f'{path_results}/CICE/CICE_XMXL_video/CICE_XMXL_{run}_no_{i:04d}')
    
    return