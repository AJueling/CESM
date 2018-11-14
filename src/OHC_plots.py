import numpy as np
import xarray as xr
import cmocean
import seaborn as sns
import matplotlib.pyplot as plt

from OHC import OHC_detrend_levels
from grid import create_tdepth, create_dz_mean
from paths import path_results
from constants import km

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
        plt.plot(ds.time/365,
                 (ds.OHC_global-ds.OHC_global[0])/1e21,
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
        plt.plot(ds.time/365,
                 (ds.OHC_global-ds.OHC_global.shift(time=1)).rolling({'time':10}).mean()/1e21,
                 c=colors[i], label=f'{labels[i]}', lw=lws[i])
    plt.text(.5,.9, f'{run.upper()}', ha='center', transform=ax.transAxes, fontsize=16)
    plt.text(.98,.02, '10 year running mean', ha='right', transform=ax.transAxes, fontsize=14)
#     plt.legend(fontsize=16, ncol=2)
    plt.xlabel('time [years]', fontsize=16)
    plt.ylabel(f'OHC($t$)-OHC($t-1$) [ZJ]', fontsize=16)
    plt.savefig(f'{path_results}/OHC/OHC_global_integrals_regional_diff_{run}')

    
def plot_global_integrals_detr(dss, run):
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
        qf = np.polyfit(ds.time , ds.OHC_global, 2)
        detr = ds.OHC_global - (qf[0]*ds.time**2 + qf[1]*ds.time + qf[2])
        plt.plot(ds.time/365,
                 (detr.rolling({'time':10}).mean())/1e21,
                 c=colors[i], label=f'{labels[i]}', lw=lws[i])
    plt.text(.5,.9, f'{run.upper()}', ha='center', transform=ax.transAxes, fontsize=16)
    plt.text(.98,.02, '10 year running mean', ha='right', transform=ax.transAxes, fontsize=14)
#     plt.legend(fontsize=16, ncol=2)
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
    times = np.arange(n)
    
    levels_trend = OHC_detrend_levels(ds.OHC_global_levels, detrend=detrend)

    X       , Y         = np.meshgrid(times, -tdepth/1e3)   # full depth
    X_detail, Y_detail  = np.meshgrid(times, -tdepth[:20])  # ca. upper 600 m
    
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
    times = np.arange(n)
    
    names = ['Southern', 'Atlantic', 'Pacific', 'Indian']
    X, Y  = np.meshgrid(times, -tdepth[:20])  # ca. upper 600 m
    
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
        


