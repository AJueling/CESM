import sys
sys.path.append("..")
import scipy as sp
import numpy as np
import xarray as xr
import seaborn as sns
import cmocean
import cartopy
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from maps import map_robinson, map_eq_earth, rect_polygon, make_map, regr_map
from paths import path_results, path_samoc
from timeseries import lowpass, chebychev
from xr_regression import xr_lintrend, lag_linregress_3D


def SST_regr_standard(index):
    """ SST regression for a given SST index
    uses 100 years of RCP run, (last) 200 years of CTRL, LPD, LPI runs
    """
    assert index in ['AMO', 'TPI', 'SOM']
    
    SST_yrly_detr_ctrl = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_detr_ctrl.nc', decode_times=False)
    SST_yrly_detr_rcp  = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_detr_rcp.nc' , decode_times=False)
    SST_yrly_detr_lpd  = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_detr_lpd.nc' , decode_times=False)[-200:,:,:]
    SST_yrly_detr_lpi  = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_detr_lpi.nc' , decode_times=False)[-200:,:,:]
    
    if index=='AMO':
        AMO_ctrl = xr.open_dataarray(f'{path_samoc}/SST/AMO_yrly_ctrl.nc', decode_times=False)
        AMO_rcp  = xr.open_dataarray(f'{path_samoc}/SST/AMO_yrly_rcp.nc' , decode_times=False)
        AMO_lpd  = xr.open_dataarray(f'{path_samoc}/SST/AMO_yrly_lpd.nc' , decode_times=False)[-200:]
        AMO_lpi  = xr.open_dataarray(f'{path_samoc}/SST/AMO_yrly_lpi.nc' , decode_times=False)[-200:]
        index_filt_detr_ctrl = chebychev(AMO_ctrl, 13) - xr_lintrend(chebychev(AMO_ctrl,13))
        index_filt_detr_rcp  = chebychev(AMO_rcp , 13) - xr_lintrend(chebychev(AMO_rcp ,13))
        index_filt_detr_lpd  = chebychev(AMO_lpd , 13) - xr_lintrend(chebychev(AMO_lpd ,13))
        index_filt_detr_lpi  = chebychev(AMO_lpi , 13) - xr_lintrend(chebychev(AMO_lpi ,13))
      
    elif index=='TPI':
        TPI_ctrl = xr.open_dataarray(f'{path_samoc}/SST/TPI_ctrl.nc', decode_times=False)
        TPI_rcp  = xr.open_dataarray(f'{path_samoc}/SST/TPI_rcp.nc' , decode_times=False)
        TPI_lpd  = xr.open_dataarray(f'{path_samoc}/SST/TPI_lpd.nc' , decode_times=False)[-200:]
        TPI_lpi  = xr.open_dataarray(f'{path_samoc}/SST/TPI_lpi.nc' , decode_times=False)[-200:]
        index_filt_detr_ctrl = chebychev(TPI_ctrl, 13) - xr_lintrend(chebychev(TPI_ctrl,13))
        index_filt_detr_rcp  = chebychev(TPI_rcp , 13) - xr_lintrend(chebychev(TPI_rcp ,13))
        index_filt_detr_lpd  = chebychev(TPI_lpd , 13) - xr_lintrend(chebychev(TPI_lpd ,13))
        index_filt_detr_lpi  = chebychev(TPI_lpi , 13) - xr_lintrend(chebychev(TPI_lpi ,13))

    elif index=='SOM':
        index_filt_detr_ctrl = xr.open_dataarray(f'{path_samoc}/SST/SOM_cheb13_index_ctrl.nc', decode_times=False)
        index_filt_detr_rcp  = xr.open_dataarray(f'{path_samoc}/SST/SOM_cheb13_index_rcp.nc' , decode_times=False)
        index_filt_detr_lpd  = xr.open_dataarray(f'{path_samoc}/SST/SOM_cheb13_index_lpd.nc' , decode_times=False)
        index_filt_detr_lpi  = xr.open_dataarray(f'{path_samoc}/SST/SOM_cheb13_index_lpi.nc' , decode_times=False)
        
#     elif index=='PDO':
        
        
    ds_ctrl = lag_linregress_3D(index_filt_detr_ctrl, SST_yrly_detr_ctrl)
    ds_rcp  = lag_linregress_3D(index_filt_detr_rcp , SST_yrly_detr_rcp )
    ds_lpd  = lag_linregress_3D(index_filt_detr_lpd , SST_yrly_detr_lpd )
    ds_lpi  = lag_linregress_3D(index_filt_detr_lpi , SST_yrly_detr_lpi )
    
    ds_ctrl.to_netcdf(f'{path_samoc}/SST/{index}_regr_ctrl.nc')
    ds_rcp .to_netcdf(f'{path_samoc}/SST/{index}_regr_rcp.nc' )
    ds_lpd .to_netcdf(f'{path_samoc}/SST/{index}_regr_lpd.nc' )
    ds_lpi .to_netcdf(f'{path_samoc}/SST/{index}_regr_lpi.nc' )
    
    return


def SST_regr_lpd(index):
    """
    as SST_regr_standard but then for first and last 200 years, as well as full 412 years of LPD run
    """
    SST_yrly_detr_lpd  = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_detr_lpd.nc' , decode_times=False)
    
    if index=='AMO':
        AMO_lpd_200_1  = xr.open_dataarray(f'{path_results}/SST/AMO_yrly_lpd.nc' , decode_times=False)[:200] 
        AMO_lpd_200_2  = xr.open_dataarray(f'{path_results}/SST/AMO_yrly_lpd.nc' , decode_times=False)[-200:]
        AMO_lpd_412    = xr.open_dataarray(f'{path_results}/SST/AMO_yrly_lpd.nc' , decode_times=False)[:]    
        index_filt_detr_lpd_200_1  = chebychev(AMO_lpd_200_1, 13) - xr_lintrend(chebychev(AMO_lpd_200_1, 13))
        index_filt_detr_lpd_200_2  = chebychev(AMO_lpd_200_2, 13) - xr_lintrend(chebychev(AMO_lpd_200_2, 13))
        index_filt_detr_lpd_412    = chebychev(AMO_lpd_412  , 13) - xr_lintrend(chebychev(AMO_lpd_412  , 13))

    elif index=='TPI':
        TPI_lpd_200_1  = xr.open_dataarray(f'{path_results}/SST/TPI_lpd.nc' , decode_times=False)[:200] 
        TPI_lpd_200_2  = xr.open_dataarray(f'{path_results}/SST/TPI_lpd.nc' , decode_times=False)[-200:]
        TPI_lpd_412    = xr.open_dataarray(f'{path_results}/SST/TPI_lpd.nc' , decode_times=False)[:]    
        index_filt_detr_lpd_200_1  = chebychev(TPI_lpd_200_1, 13) - xr_lintrend(chebychev(TPI_lpd_200_1, 13))
        index_filt_detr_lpd_200_2  = chebychev(TPI_lpd_200_2, 13) - xr_lintrend(chebychev(TPI_lpd_200_2, 13))
        index_filt_detr_lpd_412    = chebychev(TPI_lpd_412  , 13) - xr_lintrend(chebychev(TPI_lpd_412  , 13))
        
    elif index=='SOM':
        index_filt_detr_lpd_200_1  = xr.open_dataarray(f'{path_results}/SST/SOM_cheb13_index_lpd.nc' , decode_times=False)[:200] 
        index_filt_detr_lpd_200_2  = xr.open_dataarray(f'{path_results}/SST/SOM_cheb13_index_lpd.nc' , decode_times=False)[-200:]
        index_filt_detr_lpd_412    = xr.open_dataarray(f'{path_results}/SST/SOM_cheb13_index_lpd.nc' , decode_times=False)[:]    
        
    ds_lpd_200_1  = lag_linregress_3D(index_filt_detr_lpd_200_1, SST_yrly_detr_lpd[:200] )
    ds_lpd_200_2  = lag_linregress_3D(index_filt_detr_lpd_200_2, SST_yrly_detr_lpd[-200:])
    ds_lpd_412    = lag_linregress_3D(index_filt_detr_lpd_412  , SST_yrly_detr_lpd[:]    )
    
    ds_lpd_200_1.to_netcdf(f'{path_results}/SST/{index}_regr_lpd_200_1.nc')
    ds_lpd_200_2.to_netcdf(f'{path_results}/SST/{index}_regr_lpd_200_2.nc')
    ds_lpd_412  .to_netcdf(f'{path_results}/SST/{index}_regr_lpd_412.nc'  )

    
def SST_regr_lpi(index):
    """
    as SST_regr_standard, but then for first and last 800 years, as well as full 1480 years of LPI data
    """
    SST_yrly_detr_lpi  = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_detr_lpi.nc' , decode_times=False)
    
    if index=='AMO':
        AMO_lpi_800_1  = xr.open_dataarray(f'{path_results}/SST/AMO_yrly_lpi.nc' , decode_times=False)[:800] 
        AMO_lpi_800_2  = xr.open_dataarray(f'{path_results}/SST/AMO_yrly_lpi.nc' , decode_times=False)[-800:]
        AMO_lpi_1480   = xr.open_dataarray(f'{path_results}/SST/AMO_yrly_lpi.nc' , decode_times=False)[:]    
        index_filt_detr_lpi_800_1  = chebychev(AMO_lpi_800_1, 13) - xr_lintrend(chebychev(AMO_lpi_800_1, 13))
        index_filt_detr_lpi_800_2  = chebychev(AMO_lpi_800_2, 13) - xr_lintrend(chebychev(AMO_lpi_800_2, 13))
        index_filt_detr_lpi_1480   = chebychev(AMO_lpi_1480 , 13) - xr_lintrend(chebychev(AMO_lpi_1480 , 13))

    elif index=='TPI':
        TPI_lpi_800_1  = xr.open_dataarray(f'{path_results}/SST/TPI_lpi.nc' , decode_times=False)[:800] 
        TPI_lpi_800_2  = xr.open_dataarray(f'{path_results}/SST/TPI_lpi.nc' , decode_times=False)[-800:]
        TPI_lpi_1480   = xr.open_dataarray(f'{path_results}/SST/TPI_lpi.nc' , decode_times=False)[:]    
        index_filt_detr_lpi_800_1  = chebychev(TPI_lpi_800_1, 13) - xr_lintrend(chebychev(TPI_lpi_800_1, 13))
        index_filt_detr_lpi_800_2  = chebychev(TPI_lpi_800_2, 13) - xr_lintrend(chebychev(TPI_lpi_800_2, 13))
        index_filt_detr_lpi_1480   = chebychev(TPI_lpi_1480 , 13) - xr_lintrend(chebychev(TPI_lpi_1480 , 13))
    
    elif index=='SOM':
        index_filt_detr_lpi_800_1 = xr.open_dataarray(f'{path_results}/SST/SOM_cheb13_lpi.nc' , decode_times=False)[:800] 
        index_filt_detr_lpi_800_2 = xr.open_dataarray(f'{path_results}/SST/SOM_cheb13_lpi.nc' , decode_times=False)[-800:]
        index_filt_detr_lpi_1480  = xr.open_dataarray(f'{path_results}/SST/SOM_cheb13_lpi.nc' , decode_times=False)[:]
        
    ds_lpi_800_1  = lag_linregress_3D(index_filt_detr_lpi_800_1, SST_yrly_detr_lpi[:800] )
    ds_lpi_800_2  = lag_linregress_3D(index_filt_detr_lpi_800_2, SST_yrly_detr_lpi[-800:])
    ds_lpi_1480   = lag_linregress_3D(index_filt_detr_lpi_1480 , SST_yrly_detr_lpi[:]    )
    
    ds_lpi_800_1.to_netcdf(f'{path_results}/SST/{index}_regr_lpi_800_1.nc')
    ds_lpi_800_2.to_netcdf(f'{path_results}/SST/{index}_regr_lpi_800_2.nc')
    ds_lpi_1480 .to_netcdf(f'{path_results}/SST/{index}_regr_lpi_1480.nc' )


def regr_map_standard(run, index):
    assert run in ['ctrl','rcp','lpd', 'lpi']
    assert index in ['SOM', 'AMO', 'TPI']
    SST_regr_standard(index)
    ds = xr.open_dataset(f'{path_results}/SST/{index}_regr_{run}.nc', decode_times=False)
    regr_map(ds=ds, index=index, run=run)
    return


def regr_map_diff_times(run, index):
    assert run in ['lpd', 'lpi']
    assert index in ['SOM','AMO', 'TPI']
    if run=='lpd':
        times = ['200_1', '200_2', '412']
        SST_regr_lpd(index)
    elif run=='lpi':
        times = ['800_1', '800_2', '1480']
        SST_regr_lpi(index)
    for i in range(3):
        ds = xr.open_dataset(f'{path_results}/SST/{index}_regr_{run}_{times[i]}.nc', decode_times=False)
        regr_map(ds=ds, index=index, run=run, fn=times[i])
    return


if __name__=="__main__":
    for index in ['AMO', 'TPI', 'SOM']:
        print(index)
        for run in ['ctrl', 'rcp', 'lpi', 'lpd']:
            print(run)
            regr_map_standard(run, index)
            if run in ['lpd', 'lpi']:
                regr_map_diff_times(run, index)