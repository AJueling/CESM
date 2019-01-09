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


# from OHC import t2da, t2ds
# from SST import SST_index, EOF_SST_analysis
from maps import map_robinson, map_eq_earth, rect_polygon, make_map, regr_map
# from grid import find_array_idx
from paths import path_results, path_samoc#, file_ex_ocn_ctrl, file_ex_ocn_rect
# from regions import boolean_mask, SOM_area, Nino12, Nino34, global_ocean,\
#                     gl_ocean_rect, NPacific_mask_rect,\
#                     Nino12_low, Nino34_low, TexT_mask, AMO_mask, SST_index_bounds
# from plotting import shifted_color_map, discrete_cmap
from timeseries import lowpass, chebychev
from xr_regression import xr_lintrend, lag_linregress_3D


def SST_regr_standard(index):
    """
    selects last 200 years of data if there are more than 200 years
    """
    SST_yrly_detr_ctrl = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_detr_ctrl.nc', decode_times=False)
    SST_yrly_detr_rcp  = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_detr_rcp.nc' , decode_times=False)
    SST_yrly_detr_lpd  = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_detr_lpd.nc' , decode_times=False)[-200:,:,:]
    SST_yrly_detr_lpi  = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_detr_lpi.nc' , decode_times=False)[-200:,:,:]
    
    if index=='AMO':
        AMO_ctrl = xr.open_dataarray(f'{path_results}/SST/AMO_yrly_ctrl.nc', decode_times=False)
        AMO_rcp  = xr.open_dataarray(f'{path_results}/SST/AMO_yrly_rcp.nc' , decode_times=False)
        AMO_lpd  = xr.open_dataarray(f'{path_results}/SST/AMO_yrly_lpd.nc' , decode_times=False)[-200:]
        AMO_lpi  = xr.open_dataarray(f'{path_results}/SST/AMO_yrly_lpi.nc' , decode_times=False)[-200:]
        index_filt_detr_ctrl = chebychev(AMO_ctrl, 13) - xr_lintrend(chebychev(AMO_ctrl,13))
        index_filt_detr_rcp  = chebychev(AMO_rcp , 13) - xr_lintrend(chebychev(AMO_rcp ,13))
        index_filt_detr_lpd  = chebychev(AMO_lpd , 13) - xr_lintrend(chebychev(AMO_lpd ,13))
        index_filt_detr_lpi  = chebychev(AMO_lpi , 13) - xr_lintrend(chebychev(AMO_lpi ,13))
      
    elif index=='TPI':
        TPI_ctrl = xr.open_dataarray(f'{path_results}/SST/TPI_ctrl.nc', decode_times=False)
        TPI_rcp  = xr.open_dataarray(f'{path_results}/SST/TPI_rcp.nc' , decode_times=False)
        TPI_lpd  = xr.open_dataarray(f'{path_results}/SST/TPI_lpd.nc' , decode_times=False)[-200:]
        TPI_lpi  = xr.open_dataarray(f'{path_results}/SST/TPI_lpi.nc' , decode_times=False)[-200:]
        index_filt_detr_ctrl = chebychev(TPI_ctrl, 13) - xr_lintrend(chebychev(TPI_ctrl,13))
        index_filt_detr_rcp  = chebychev(TPI_rcp , 13) - xr_lintrend(chebychev(TPI_rcp ,13))
        index_filt_detr_lpd  = chebychev(TPI_lpd , 13) - xr_lintrend(chebychev(TPI_lpd ,13))
        index_filt_detr_lpi  = chebychev(TPI_lpi , 13) - xr_lintrend(chebychev(TPI_lpi ,13))

    elif index=='SOM':
        SOM_ctrl = xr.open_dataarray(f'{path_results}/SST/SOM_index_ctrl.nc', decode_times=False)
        SOM_rcp  = xr.open_dataarray(f'{path_results}/SST/SOM_index_rcp.nc' , decode_times=False)
        SOM_lpd  = xr.open_dataarray(f'{path_results}/SST/SOM_index_lpd.nc' , decode_times=False)
        SOM_lpi  = xr.open_dataarray(f'{path_results}/SST/SOM_index_lpi.nc' , decode_times=False)
        index_filt_detr_ctrl = chebychev(SOM_ctrl, 13) - xr_lintrend(chebychev(SOM_ctrl,13))
        index_filt_detr_rcp  = chebychev(SOM_rcp , 13) - xr_lintrend(chebychev(SOM_rcp ,13))
        index_filt_detr_lpd  = chebychev(SOM_lpd , 13) - xr_lintrend(chebychev(SOM_lpd ,13))
        index_filt_detr_lpi  = chebychev(SOM_lpi , 13) - xr_lintrend(chebychev(SOM_lpi ,13))
        
    ds_ctrl = lag_linregress_3D(index_filt_detr_ctrl, SST_yrly_detr_ctrl)
    ds_rcp  = lag_linregress_3D(index_filt_detr_rcp , SST_yrly_detr_rcp )
    ds_lpd  = lag_linregress_3D(index_filt_detr_lpd , SST_yrly_detr_lpd )
    ds_lpi  = lag_linregress_3D(index_filt_detr_lpi , SST_yrly_detr_lpi )
    
    ds_ctrl.to_netcdf(f'{path_results}/SST/{index}_regr_ctrl.nc')
    ds_rcp .to_netcdf(f'{path_results}/SST/{index}_regr_rcp.nc' )
    ds_lpd .to_netcdf(f'{path_results}/SST/{index}_regr_lpd.nc' )
    ds_lpi .to_netcdf(f'{path_results}/SST/{index}_regr_lpi.nc' )
    return


def SST_regr_lpd(index):
    """
    selects last 200 years of data if there are more than 200 years
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
        SOM_lpd_200_1  = xr.open_dataarray(f'{path_results}/SST/SOM_index_lpd.nc' , decode_times=False)[:200] 
        SOM_lpd_200_2  = xr.open_dataarray(f'{path_results}/SST/SOM_index_lpd.nc' , decode_times=False)[-200:]
        SOM_lpd_412    = xr.open_dataarray(f'{path_results}/SST/SOM_index_lpd.nc' , decode_times=False)[:]    
        index_filt_detr_lpd_200_1  = chebychev(SOM_lpd_200_1, 13) - xr_lintrend(chebychev(SOM_lpd_200_1, 13))
        index_filt_detr_lpd_200_2  = chebychev(SOM_lpd_200_2, 13) - xr_lintrend(chebychev(SOM_lpd_200_2, 13))
        index_filt_detr_lpd_412    = chebychev(SOM_lpd_412  , 13) - xr_lintrend(chebychev(SOM_lpd_412  , 13))
        
        
    ds_lpd_200_1  = lag_linregress_3D(index_filt_detr_lpd_200_1, SST_yrly_detr_lpd[:200] )
    ds_lpd_200_2  = lag_linregress_3D(index_filt_detr_lpd_200_2, SST_yrly_detr_lpd[-200:])
    ds_lpd_412    = lag_linregress_3D(index_filt_detr_lpd_412  , SST_yrly_detr_lpd[:]    )
    
    ds_lpd_200_1.to_netcdf(f'{path_results}/SST/{index}_regr_lpd_200_1.nc')
    ds_lpd_200_2.to_netcdf(f'{path_results}/SST/{index}_regr_lpd_200_2.nc')
    ds_lpd_412  .to_netcdf(f'{path_results}/SST/{index}_regr_lpd_412.nc'  )

    
def SST_regr_lpi(index):
    """
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
        SOM_lpi_800_1  = xr.open_dataarray(f'{path_results}/SST/SOM_lpi.nc' , decode_times=False)[:800] 
        SOM_lpi_800_2  = xr.open_dataarray(f'{path_results}/SST/SOM_lpi.nc' , decode_times=False)[-800:]
        SOM_lpi_1480   = xr.open_dataarray(f'{path_results}/SST/SOM_lpi.nc' , decode_times=False)[:]    
        index_filt_detr_lpi_800_1  = chebychev(SOM_lpi_800_1, 13) - xr_lintrend(chebychev(SOM_lpi_800_1, 13))
        index_filt_detr_lpi_800_2  = chebychev(SOM_lpi_800_2, 13) - xr_lintrend(chebychev(SOM_lpi_800_2, 13))
        index_filt_detr_lpi_1480   = chebychev(SOM_lpi_1480 , 13) - xr_lintrend(chebychev(SOM_lpi_1480 , 13))
        
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