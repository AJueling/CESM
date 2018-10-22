import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt


def map_robinson(xr_DataArray, cmap, minv, maxv, label, filename=None):
    fig = plt.figure(figsize=(8,5))
    ax  = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_position([.02,.05,.96,.93])
    cax, kw = mpl.colorbar.make_axes(ax,location='bottom',pad=0.03,shrink=0.8)
    im = ax.pcolormesh(xr_DataArray.lon,
                       xr_DataArray.lat,
                       xr_DataArray.values,
                       cmap=cmap,
                       vmin=minv, vmax=maxv,
                       transform=ccrs.PlateCarree() )
    ax.coastlines()
    cbar = fig.colorbar(im, cax=cax, extend='both', **kw)
    cbar.ax.tick_params(labelsize=14)
    label = cbar.set_label(label, size=16)
    if filename!=None: plt.savefig(filename)
    return fig


def map_ocn_robinson(xr_DataArray, cmap, minv, maxv, label, filename=None, grid='T'):
    fig = plt.figure(figsize=(8,5))
    ax  = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_position([.02,.05,.96,.93])
    cax, kw = mpl.colorbar.make_axes(ax,location='bottom',pad=0.03,shrink=0.8)
    if grid=='T':
        im = ax.pcolormesh(xr_DataArray.TLONG,
                           xr_DataArray.TLAT,
                           xr_DataArray.values,
                           cmap=cmap,
                           vmin=minv, vmax=maxv,
                           transform=ccrs.PlateCarree() )
        ax.coastlines()
    elif grid=='U':
        im = ax.pcolormesh(xr_DataArray.ULONG,
                           xr_DataArray.ULAT,
                           xr_DataArray.values,
                           cmap=cmap,
                           vmin=minv, vmax=maxv,
                           transform=ccrs.PlateCarree() )
        ax.add_feature(cartopy.feature.LAND, zorder=2, edgecolor='black', facecolor='w')
    cbar = fig.colorbar(im, cax=cax, extend='both', **kw)
    cbar.ax.tick_params(labelsize=14)
    label = cbar.set_label(label, size=16)
    if filename!=None: plt.savefig(filename)
    return fig