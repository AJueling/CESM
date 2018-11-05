import numpy as np
import xarray as xr
import seaborn as sns
import matplotlib.pyplot as plt

from grid import create_tdepth
from paths import path_results

tdepth = create_tdepth(domain='ocn')
colors = ['k', 'C0', 'C1', 'C2', 'C3', 'C4']
labels = ['Global', 'Atlantic', 'Pacific', 'Indian', 'Southern', 'Mediterranean']


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
        plt.plot((ds.OHC_global-ds.OHC_global[0])/1e21, c=colors[i], lw=2, label=labels[i])
    plt.text(.9,.7, f'{run.upper()}', transform=ax.transAxes, fontsize=16)
    plt.xlabel('time [years]', fontsize=16)
    plt.ylabel('OHC [ZJ]', fontsize=16)
    plt.legend(fontsize=16)
    plt.savefig(f'{path_results}/OHC/OHC_global_integrals_regional_{run}')
    

    
def plot_global_integrals_detrended(dss, run):
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
        plt.plot((detr.rolling({'time': 5}).mean())/1e21, c=colors[i], label=f'{labels[i]}', lw=2)
    plt.text(.9,.9, f'{run.upper()}', transform=ax.transAxes, fontsize=16)
    plt.legend(fontsize=16, ncol=2)
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
    assert run in ['rcp', 'ctrl']
    
    f = plt.figure(figsize=(8,5))
    ax = f.add_axes([0.13,0.13,.85,.85])
    plt.tick_params(labelsize=14)
    plt.axvline(0, c='k', lw=.5)
    for i in range(len(das)):
        plt.plot(das[i], -tdepth/1e3, c=colors[i] , lw=2, label=f'{labels[i]}')
    plt.text(.05, .9, f'{run.upper()}', transform=ax.transAxes, fontsize=16)
    plt.xlim((-1e18,4e18))
    plt.xlabel('OHC trend [J/m/yr]', fontsize=16)
    plt.ylabel('depth [km]', fontsize=16)
    plt.legend(fontsize=16)
    plt.tight_layout()
    plt.savefig(f'{path_results}/OHC/OHC_global_levels_regional_{run}')
    

    
def Hovmoeller_global_depth(ds):
    """ plots full depth / surface detail Hovmoeller diagrams of OHC(z) """
    
    