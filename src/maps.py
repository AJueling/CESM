import numpy as np
import cartopy
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt


def map_robinson(xa, domain, cmap, minv, maxv, label, filename=None):
    
    assert domain in ['atm', 'ocn_T', 'ocn_U']
    
    fig = plt.figure(figsize=(8,5))
    ax  = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_position([.02,.05,.96,.93])
    cax, kw = mpl.colorbar.make_axes(ax,location='bottom',pad=0.03,shrink=0.8)
    
    if domain=='atm':
        im = ax.pcolormesh(xa.lon,
                           xa.lat,
                           xa.values,
                           cmap=cmap,
                           vmin=minv, vmax=maxv,
                           transform=ccrs.PlateCarree() )
        ax.coastlines()
        
    elif domain=='ocn_T':
        im = ax.pcolormesh(xa.TLONG,
                           xa.TLAT,
                           xa.values,
                           cmap=cmap,
                           vmin=minv, vmax=maxv,
                           transform=ccrs.PlateCarree() )
        ax.add_feature(cartopy.feature.LAND, zorder=2, edgecolor='black', facecolor='w')
        
    elif domain=='ocn_U':
        im = ax.pcolormesh(xa.ULONG,
                           xa.ULAT,
                           xa.values,
                           cmap=cmap,
                           vmin=minv, vmax=maxv,
                           transform=ccrs.PlateCarree() )
        ax.add_feature(cartopy.feature.LAND, zorder=2, edgecolor='black', facecolor='w')
    
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180])
    
    cbar = fig.colorbar(im, cax=cax, extend='both', **kw)
    cbar.ax.tick_params(labelsize=14)
    label = cbar.set_label(label, size=16)
    if filename!=None: plt.savefig(filename)
    return fig, ax


def map_ocn_robinson(xr_DataArray, cmap, minv, maxv, label, filename=None, grid='T'):
    fig = plt.figure(figsize=(8,5))
    ax  = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    ax.set_position([.02,.05,.96,.93])
    cax, kw = mpl.colorbar.make_axes(ax,location='bottom',pad=0.03,shrink=0.8)
    
    cbar = fig.colorbar(im, cax=cax, extend='both', **kw)
    cbar.ax.tick_params(labelsize=14)
    label = cbar.set_label(label, size=16)
    if filename!=None: plt.savefig(filename)
    return fig