import numpy as np
import xarray as xr
import cmocean
import seaborn as sns
import cmocean
import matplotlib as mpl
import matplotlib.pyplot as plt

from OHC import OHC_detrend_levels
from maps import map_eq_earth
from grid import create_tdepth, create_dz_mean
from paths import path_results, path_samoc
from plotting import discrete_cmap, shifted_color_map
from constants import km, spy
from timeseries import lowpass
from xr_regression import ocn_field_regression, xr_quadtrend

tdepth = create_tdepth(domain='ocn')
colors = ['k', 'C0', 'C1', 'C2', 'C3', 'C4']
labels = ['Global', 'Atlantic', 'Pacific', 'Indian', 'Southern', 'Mediterranean']
lws    = [2.5,1.5,1.5,1,1.5,1]


def plot_global_integrals(dss, run):
    """"""
    n = len(dss)
    assert n<=6
    assert run in ['rcp', 'ctrl']
    
    f = plt.figure(figsize=(8,5))
    ax = f.add_axes([0.13,0.13,.85,.85])
    plt.tick_params(labelsize=14)
    plt.axhline(0, c='k', lw=.5)
    for i, ds in enumerate(dss):
        plt.plot(ds.time/365, (ds.OHC_global-ds.OHC_global[0])/1e21,
                 c=colors[i], lw=lws[i], label=labels[i])
    plt.text(.5,.9, f'{run.upper()}', ha='center', transform=ax.transAxes, fontsize=16)
    plt.xlabel('time [years]', fontsize=16)
    plt.ylabel('OHC [ZJ]', fontsize=16)
    plt.legend(fontsize=16)
    plt.savefig(f'{path_results}/OHC/OHC_global_integrals_regional_{run}')
    
    
def plot_global_integrals_diff(dss, run):
    """
    
    input:
    dss .. list of datasets
    """
    n = len(dss)
    assert n<=6
    assert run in ['rcp', 'ctrl']
    
    f = plt.figure(figsize=(8,5))
    ax = f.add_axes([0.13,0.13,.85,.85])
    plt.tick_params(labelsize=14)
    plt.axhline(0, c='k', lw=.5)
    for i, ds in enumerate(dss):
#         plt.plot((ds.OHC_global-ds.OHC_global.shift(time=1))/1e21, c=colors[i], lw=.5)
        plt.plot(ds.time[1:]/365, lowpass((ds.OHC_global-ds.OHC_global.shift(time=1))[1:], 10)/1e21,
                 c=colors[i], label=f'{labels[i]}', lw=lws[i])
    plt.text(.5,.9, f'{run.upper()}', ha='center', transform=ax.transAxes, fontsize=16)
    plt.text(.98,.02, '10 lowpass filtered', ha='right', transform=ax.transAxes, fontsize=14)
#     plt.legend(fontsize=16, ncol=2)
    plt.xlabel('time [years]', fontsize=16)
    plt.ylabel(f'OHC($t$)-OHC($t-1$) [ZJ]', fontsize=16)
    plt.savefig(f'{path_results}/OHC/OHC_global_integrals_regional_diff_{run}')

    
def plot_global_integrals_detr(dss, run):
    """ quadratically detrended
    
    input:
    dss .. list of datasets
    """
    n = len(dss)
    assert n<=6
    assert run in ['rcp', 'ctrl']
    
    f = plt.figure(figsize=(8,5))
    ax = f.add_axes([0.13,0.13,.85,.85])
    plt.tick_params(labelsize=14)
    plt.axhline(0, c='k', lw=.5)
    for i, ds in enumerate(dss):
        qf = np.polyfit(ds.time , ds.OHC_global, 2)
        detr = ds.OHC_global - xr_quadtrend(ds.OHC_global)
        plt.plot(ds.time/365, lowpass(detr, 10)/1e21,
                 c=colors[i], label=f'{labels[i]}', lw=lws[i])
    plt.text(.5,.9, f'{run.upper()}', ha='center', transform=ax.transAxes, fontsize=16)
    plt.text(.98,.02, '10 lowpass filtered', ha='right', transform=ax.transAxes, fontsize=14)
    plt.xlabel('time [years]', fontsize=16)
    plt.ylabel('quad. detrended OHC [ZJ]', fontsize=16)
    plt.savefig(f'{path_results}/OHC/OHC_global_integrals_regional_detr_{run}')

    
    

def plot_levels_trend(das, run):
    """ plots growth rates per level
    
    input:
    das .. list of xr DataArrays
    run .. rcp or ctrl
    """
    n = len(das)
    assert n<=6
    for i in range(n):  assert len(das[i])==len(tdepth)
    assert run in ['rcp', 'ctrl', 'rcp-ctrl']
    
    f = plt.figure(figsize=(8,5))
    ax = f.add_axes([0.10,0.13,.88,.85])
    plt.tick_params(labelsize=14)
    plt.axvline(0, c='k', lw=.5)
    for i in range(len(das)):
        plt.plot(das[i], -tdepth/1e3, c=colors[i] , lw=2, label=f'{labels[i]}')
    plt.text(.02, .9, f'{run.upper()}', transform=ax.transAxes, fontsize=16)
    plt.xlim((-1e18,4e18))
    plt.xlabel('OHC trend [J/m/yr]', fontsize=16)
    plt.ylabel('depth [km]', fontsize=16)
    plt.legend(fontsize=16)
    plt.tight_layout()
    plt.savefig(f'{path_results}/OHC/OHC_global_levels_regional_{run}')
    

    
def Hovmoeller_global_depth(ds, detrend='lin', fn=None):
    """ plots full depth / surface detail Hovmoeller diagrams of OHC(z) """
    assert detrend in ['lin', 'quad']
    assert 'OHC_global_levels' in ds

    dz_mean = create_dz_mean(domain='ocn')

    n = len(ds.OHC_global_levels.time)
    
    levels_trend = OHC_detrend_levels(ds.OHC_global_levels, detrend=detrend)

    X       , Y         = np.meshgrid(ds.time/365, -tdepth/1e3)   # full depth
    X_detail, Y_detail  = np.meshgrid(ds.time/365, -tdepth[:20])  # ca. upper 600 m
    
    f, ax = plt.subplots(2, 1, figsize=(8,5), sharex=True)

    maxv = .2
    ax[0].pcolormesh(X_detail, Y_detail, levels_trend[:,:20].T/1e21, cmap=cmocean.cm.balance, vmin=-maxv, vmax=maxv)
    ax[1].pcolormesh(X       , Y       , levels_trend.T/1e21       , cmap=cmocean.cm.balance, vmin=-maxv, vmax=maxv)
    ax[1].set_xlabel('time [years]', fontsize=16)
    for i in range(2):
        ax[i].tick_params(labelsize=14)
    ax[0].set_yticks([0,-200,-400,-600])
    ax[1].set_yticks([0,-2,-4,-6])
    ax[0].set_ylabel('depth [m]' , fontsize=16)
    ax[1].set_ylabel('depth [km]', fontsize=16)
    ax[1].text(.01, .05, f'{detrend}. detr.', ha='left', transform=ax[1].transAxes, fontsize=14, color='grey')
    plt.tight_layout()
    f.align_labels()
    if fn!=None:  plt.savefig(fn)

        
        
def Hovmoeller_basins_depth(dss, detrend='lin', fn=None):
    """  """
    assert detrend in ['lin', 'quad']
    assert len(dss)==4
    
    n = len(dss[0].OHC_global_levels.time)
    
    names = ['Southern', 'Atlantic', 'Pacific', 'Indian']
    X, Y  = np.meshgrid(dss[0].time/365, -tdepth[:20])  # ca. upper 600 m
    
    f, ax = plt.subplots(4, 1, figsize=(8,8), sharex=True)
    maxv=.1
    
    for j, ds in enumerate(dss):

        levels_trend = OHC_detrend_levels(ds.OHC_global_levels, detrend=detrend)
        
        ax[j].pcolormesh(X, Y, levels_trend[:,:20].T/1e21, cmap=cmocean.cm.balance, vmin=-maxv, vmax=maxv)

        ax[j].set_yticks([0,-200,-400,-600])
        ax[j].set_ylabel('depth [m]' , fontsize=16)
        ax[j].tick_params(labelsize=14)
        ax[j].text(.98, .05, f'{names[j]} Ocean', ha='right', transform=ax[j].transAxes, fontsize=16)
        ax[j].text(.01, .05, f'{detrend}. detr.', ha='left', transform=ax[j].transAxes, fontsize=14, color='grey')
    ax[3].set_xlabel('time [years]', fontsize=16)
    plt.tight_layout()

    
    if fn!=None:  plt.savefig(fn)
        

def OHC_vert_video(run):
    domain = 'ocn_T'

    if run=='ctrl':
        first_year = 100
        nt = 180
        ds = xr.open_dataset(f'{path_samoc}/OHC/OHC_integrals_Global_Ocean_ctrl.nc'  , decode_times=False)
        # something is still wrong in CTRL year 205
#         for ds in [ctrl]:#, ctrl_A, ctrl_P, ctrl_I, ctrl_S]:
        for field in ['OHC_global', 'OHC_global_levels', 'OHC_zonal', 'OHC_zonal_levels']:
            ds[field][105] = (ds[field].sel({'time':204*365}) +
                                           ds[field].sel({'time':206*365}) )/2
        OHC_vt   = xr.open_dataarray(f'{path_samoc}/OHC/OHC_trend_vert_int_c2.nc'    , decode_times=False)
        OHC_vt_a = xr.open_dataarray(f'{path_samoc}/OHC/OHC_trend_vert_int_a_c2.nc'  , decode_times=False)
        OHC_vt_b = xr.open_dataarray(f'{path_samoc}/OHC/OHC_trend_vert_int_b_c2.nc'  , decode_times=False)

    elif run=='rcp':
        first_year = 200
        nt = 80
        ds = xr.open_dataset(f'{path_samoc}/OHC/OHC_integrals_Global_Ocean_rcp.nc'   , decode_times=False)
        OHC_vt   = xr.open_dataarray(f'{path_samoc}/OHC/OHC_trend_vert_int_rcp.nc'   , decode_times=False)
        OHC_vt_a = xr.open_dataarray(f'{path_samoc}/OHC/OHC_trend_vert_int_a_rcp.nc' , decode_times=False)
        OHC_vt_b = xr.open_dataarray(f'{path_samoc}/OHC/OHC_trend_vert_int_b_rcp.nc' , decode_times=False)

    vert_trends = [OHC_vt, OHC_vt_a, OHC_vt_b]
    names1 = ['full_depth', 'above_100m', 'below_100m']
    names2 = ['full depth', 'above 100m', 'below 100m']
    names3 = [f'full\ndepth', f'above\n100m', f'below\n100m']

    for t in range(nt-10):
        dss = [ds.OHC_vertical[t:t+10,:,:],
               ds.OHC_vertical_above_100m[t:t+10,:,:],
               ds.OHC_vertical_below_100m[t:t+10,:,:]]

        for i, xa in enumerate(dss):
            text2  = names3[i]
            
            text1 = f'{t+first_year}-{t+first_year+10}'
            label = r'10 year OHC trend [W/m$^2$]'
            if i in [0,2]:
                minv  = -3
                maxv  = 9
                nc    = 12
            else:
                minv  = -2
                maxv  = 6
                nc    = 8
            cmap  = discrete_cmap(nc, shifted_color_map(cmocean.cm.balance, 
                                                        start=.33, midpoint=0.5,
                                                        stop=1., name='shrunk'))
            xa1 = ocn_field_regression(xa)/spy
            fn = f'{path_results}/OHC/OHC-video/OHC_trend_vert_int_{names1[i]}_map_{run}_{t}'
            f, ax = map_eq_earth(xa=xa1, domain=domain, cmap=cmap, minv=minv, maxv=maxv,
                                 label=label, filename=fn, text1=text1, text2=text2)

            # anomaly from long term trend
            text1 = f'({t+first_year}-{t+first_year+10})-\n<RCP>'
            label = r'10 year OHC trend anomaly [W/m$^2$]'
            if i in [0,2]:  # full depth, below 100 m
                minv   = -10
                maxv   = 10
                nc = 20
            else:           # above 100 m
                minv   = -5
                maxv   = 5
                nc = 10
            cmap   = discrete_cmap(nc, cmocean.cm.curl)
            xa2 = xa1 - vert_trends[i]
            fn = f'{path_results}/OHC/OHC-video/OHC_trend_vert_int_{names1[i]}_anom_map_{run}_{t}'
            f, ax = map_eq_earth(xa=xa2, domain=domain, cmap=cmap, minv=minv, maxv=maxv,
                                 label=label, filename=fn, text1=text1, text2=text2)
            
            plt.close('all')
            
            