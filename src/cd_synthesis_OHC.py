import os
import sys
sys.path.append("..")
import numpy as np
import xarray as xr
import cmocean
import matplotlib
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs


from tqdm import tqdm
from paths import path_samoc, path_results
from regions import boolean_mask

matplotlib.rc_file('../rc_file')
matplotlib.use('Agg')


class SynthesizeOHC(object):
    """ functions to generate netcdf files derived from CESM output / obs. """
    
    def __init__(self):
        return
        
        
    def plot_OHC_anomaly(self):
        """"""
        ctrl_qd = xr.open_dataset(f'{path_samoc}/OHC/OHC_integrals_ctrl_qd.nc', decode_times=False)
        lpd_qd  = xr.open_dataset(f'{path_samoc}/OHC/OHC_integrals_lpd_qd.nc' , decode_times=False)
        
        maxv = []
        for j, depths in enumerate([(0,6000), (0,100), (0,700), (700,2000)]):
            key = f'OHC_vertical_{depths[0]}_{depths[1]}m'
            maxv.append(np.max([np.abs(ctrl_qd[key]).max(), np.abs(lpd_qd[key]).max()])/4)
        print(maxv)
        
        for y in range(250):
            f, ax = plt.subplots(4, 3 , figsize=(10,10),
                                 gridspec_kw={"width_ratios":[1,1, 0.05]}, 
                                 subplot_kw=dict(projection=ccrs.EqualEarth(central_longitude=300)))
            for i, ds in enumerate([ctrl_qd, lpd_qd]):
                name = ['CTRL', 'LPD'][i]
                MASK = boolean_mask(['ocn_rect', 'ocn_low'][i], mask_nr=0)
                if i==0:   X, Y = np.meshgrid(ds.t_lon, ds.t_lat)
                else:      X, Y = ds.TLONG, ds.TLAT
                for j, depths in enumerate([(0,6000), (0,100), (0,700), (700,2000)]):
                    key = f'OHC_vertical_{depths[0]}_{depths[1]}m'
                    im = ax[j,i].pcolormesh(X, Y, ds[key][y,:,:].where(MASK),
                                            transform=ccrs.PlateCarree(),
                                            vmin=-maxv[j], vmax=maxv[j],
                                            cmap=cmocean.cm.balance)
                    ax[j,i].add_feature(cartopy.feature.LAND,
                                        zorder=2, edgecolor='black', facecolor='w')
                    if j==0:
                        year = f'{ds.time.values[y]:3.0f}'
                        ax[0,i].text(.5, 1.1, f'{name} (year {year})',
                                     transform=ax[0,i].transAxes, ha='center')
                    ax[j,i].text(.5, 1.02, ['full depth (0-6000m)', 'surface (0-100m)',
                                            'upper ocean (0-700m)', 'lower ocean (700-2000m)'][j],
                                 transform=ax[j,i].transAxes, ha='center')
                    if i==1:
                        cb = f.colorbar(im, ax=ax[j,2], orientation='vertical', label=r'OHC [J/m$^{2}$]')#, ticks=np.arange(-3e16,4e16,1e16))
                        cb.outline.set_visible(False)
            plt.savefig(f'{path_results}/OHC/OHC-video/OHC_vert_qd_ctrl_lpd_{y:03d}')
        return