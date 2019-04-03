import numpy as np
import xarray as xr
import cmocean
import cartopy
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from paths import path_results
from regions import boolean_mask, SST_index_bounds
from plotting import discrete_cmap


def map_robinson(xa, domain, cmap, minv, maxv, label, filename=None, text1=None, text2=None, rects=None, sig=None, clon=0):
    fig, ax = make_map(xa=xa, domain=domain, proj='rob', cmap=cmap, minv=minv, maxv=maxv, label=label,
                       filename=filename, text1=text1, text2=text2, rects=rects, sig=sig, clon=clon)
    return fig, ax


def map_eq_earth(xa, domain, cmap, minv, maxv, label, filename=None, text1=None, text2=None, rects=None, sig=None, clon=0):
    fig, ax = make_map(xa=xa, domain=domain, proj='ee', cmap=cmap, minv=minv, maxv=maxv, label=label,
                       filename=filename, text1=text1, text2=text2, rects=rects, sig=sig, clon=clon)
    return fig, ax


def make_map(xa, domain, proj, cmap, minv, maxv, label, filename=None, text1=None, text2=None, rects=None, sig=None, clon=0):
    """ global map (Robinson or Equal Earth projection) of xa 
    optional: significance shading, polygons, text, central longitude, file output 
    """
    assert type(xa)==xr.core.dataarray.DataArray
    assert domain in ['atm', 'ocn_T', 'ocn_U', 'ocn_rect', 'ocn_low', 'ocn_had']
    assert proj in ['ee', 'rob']
    
    fig = plt.figure(figsize=(8,5))
    if proj=='ee':
        ax  = fig.add_subplot(1, 1, 1, 
                              projection=ccrs.EqualEarth(central_longitude=clon))
    elif proj=='rob':
        ax  = fig.add_subplot(1, 1, 1,
                              projection=ccrs.Robinson(central_longitude=clon))
    ax.set_position([.02,.05,.96,.93])
    cax, kw = mpl.colorbar.make_axes(ax,location='bottom',pad=0.03,shrink=0.8)
    
    if domain in ['atm']:#, 'ocn_low']:
        lats = xa.lat
        lons = xa.lon
    elif domain in ['ocn_T', 'ocn_low']:
        lats = xa.TLAT
        lons = xa.TLONG
    elif domain=='ocn_U':
        lats = xa.ULAT
        lons = xa.ULONG
    elif domain=='ocn_rect':
        lats = xa.t_lat
        lons = xa.t_lon
    elif domain=='ocn_had':
        lats = xa.latitude
        lons = xa.longitude
        lons, lats = np.meshgrid(lons, lats)
    
    im = ax.pcolormesh(lons, lats, xa.values,
                       cmap=cmap, vmin=minv, vmax=maxv,
                       transform=ccrs.PlateCarree(),
                      )

    # significance outline
    if type(sig)==xr.core.dataarray.DataArray:
        if domain=='ocn_had':  flons, flats = lons.flatten(), lats.flatten()
        else:  flons, flats = lons.values.flatten(), lats.values.flatten()
        ax.tricontour(flons, flats, sig.values.flatten(), levels=[.5],
                      linestyles='dashed', linewidths=1.5, #cmap='gray',
                      transform=ccrs.PlateCarree(),
                     )
    
    # coastlines/land
    if domain=='atm':
        ax.coastlines()
    elif domain in ['ocn_T', 'ocn_U', 'ocn_had', 'ocn_low']:
        ax.add_feature(cartopy.feature.LAND,
                       zorder=2, edgecolor='black', facecolor='w')
        
    # text
    if text1!=None:
        ax.text(0, 1, text1, ha='left' , va='top',
                transform=ax.transAxes, fontsize=16)
    if text2!=None:
        ax.text(1, 1, text2, ha='right', va='top',
                transform=ax.transAxes, fontsize=16)
    
    # SST index polygons
    if type(rects)==np.ndarray:  rects = [rects]
    assert type(rects)==list
    for rect in rects:
        assert type(rect)==np.ndarray
        ax.add_patch(mpatches.Polygon(xy=rect,
                                      facecolor='none', edgecolor='k',
                                      linewidth=2, zorder=2,
                                      transform=ccrs.PlateCarree(), ), )
            
    # grid
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False)
    gl.ylocator = mticker.FixedLocator([-90, -60, -30, 0, 30, 60, 90])
    gl.xlocator = mticker.FixedLocator([-180, -120, -60, 0, 60, 120, 180])
    
    # colorbar
    cbar = fig.colorbar(im, cax=cax, extend='both', **kw)
    cbar.ax.tick_params(labelsize=14)
    label = cbar.set_label(label, size=16)
    
    # output
    if filename!=None: plt.savefig(filename, dpi=100)
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


def rect_polygon(extent):
    assert type(extent)==tuple
    (lonmin,lonmax,latmin,latmax) = extent
    n=50
    xs = [np.linspace(lonmin,lonmax,n), np.linspace(lonmax,lonmax,n),
          np.linspace(lonmax,lonmin,n), np.linspace(lonmin,lonmin,n)]
    ys = [np.linspace(latmin,latmin,n), np.linspace(latmin,latmax,n),
          np.linspace(latmax,latmax,n), np.linspace(latmax,latmin,n)]
    xs = [item for sublist in xs for item in sublist]
    ys = [item for sublist in ys for item in sublist]
    poly_coords = np.swapaxes(np.array([xs, ys]),0,1)
    return poly_coords


def regr_map(ds, index, run, fn=None):
    """ map of regression slope with 95% significance countours and SST index polygons """
    if run in ['ctrl', 'rcp']:   domain = 'ocn'
    elif run in ['lpd', 'lpi']:  domain = 'ocn_low'
    elif run=='had':             domain = 'ocn_had'
    MASK = boolean_mask(domain=domain, mask_nr=0)
    
    xa = ds.slope.where(MASK)
    if domain in ['ocn', 'ocn_low']:
        xa = xa.assign_coords(TLONG=ds.TLONG)
    
    if index in ['AMO', 'SOM']:
        rects = rect_polygon(SST_index_bounds(index))
        clon = 300
        nv = .4
    elif index in ['PDO', 'IPO']:
        rects = rect_polygon(SST_index_bounds(index))
        clon = 200
        nv = .4
    elif index=='TPI':
        rects = [rect_polygon(SST_index_bounds('TPI1')),
                 rect_polygon(SST_index_bounds('TPI2')),
                 rect_polygon(SST_index_bounds('TPI3')),
                ]
        clon = 200
        nv = .4
    
    print(run, np.max(ds.slope.values))
    
    # choose two-tailed 95% significance level
    # as boolean map
    sig = ds.pval#.where(MASK)
#     tail1 = np.where(sig<0.025, 1, 0)
    tail1 = np.where(sig<0.005, 1, 0)
#     tail2 = np.where(sig>0.975, 1, 0)
    tail2 = np.where(sig>99.5, 1, 0)
    sig.values = tail1 + tail2
#     if run in ['ctrl', 'rcp', 'had']:   sig = sig.where(MASK)
    
    proj = 'rob'
    cm = discrete_cmap(16, cmocean.cm.balance)    
    label ='regression slope [K/K]'
    text1 = f'SST({index})\nregr.'
    text2 = f'{run.upper()}\n{ds.first_year}-\n{ds.last_year}'
    if fn==None:  filename = f'{path_results}/SST/{index}_regr_map_{run}.png'
    else:         filename = f'{path_results}/SST/{index}_regr_map_{run}_{fn}.png'
    if run in ['ctrl', 'rcp']:   domain = 'ocn_T'
        
    f, ax = make_map(xa=xa, domain=domain, proj=proj, cmap=cm, minv=-nv, maxv=nv,
                     label=label, filename=filename, text1=text1, text2=text2,
                     rects=rects, sig=sig, clon=clon)