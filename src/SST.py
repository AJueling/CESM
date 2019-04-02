import os
import sys
import numpy as np
import xarray as xr
import datetime

from eofs.xarray import Eof

import xarray as xr

from paths import path_samoc, path_results, path_data
from regions import boolean_mask, TPI_masks, mask_box_in_region, bll_AMO, bll_SOM, bll_TPI1, bll_TPI2, bll_TPI3
from timeseries import IterateOutputCESM, deseasonalize
from xr_DataArrays import xr_AREA
from xr_regression import xr_lintrend, xr_quadtrend, lag_linregress_3D


def SST_area_average(xa_SST, AREA, AREA_index, MASK, dims=('nlat', 'nlon'), index_loc=None):
    """ calculates the average SST over an area, possibly as a time series """
    assert type(xa_SST)==xr.core.dataarray.DataArray
    assert type(AREA)==xr.core.dataarray.DataArray
    if type(index_loc)==dict:
        index = (xa_SST*AREA).where(MASK).sel(index_loc).sum(dim=dims)/AREA_index
    elif index_loc==None:
        index = (xa_SST*AREA).where(MASK).sum(dim=dims)/AREA_index
    else:
        print('kwarg `index_loc` is not given properly.')
    return index


def SST_index_from_monthly(run, index_loc, MASK):
    """ loads monthly SST data, calculated SST_index, returns raw timeseries"""
    assert run in ['ctrl', 'rcp', 'lpd', 'lpi']
    if run in ['ctrl', 'rcp']:
        domain = 'ocn'
    elif run in ['lpd', 'lpi']:
        domain = 'ocn_low'
    AREA = xr_AREA(domain)
    AREA_index = AREA.where(MASK).sum()
    
    for  i, (y, m, s) in enumerate(IterateOutputCESM(domain=domain, run=run, tavg='monthly')):
        if m==1: print(y)
        xa_SST = xr.open_dataset(s, decode_times=False).TEMP[0,0,:,:]
        SSTi = SST_index(xa_SST, AREA, index_loc, AREA_index, MASK)
        if i==0:
            new_SSTi = SSTi.copy()
        else:
            new_SSTi = xr.concat([new_SSTi,SSTi], dim='time')
            
    return new_SSTi


def EOF_SST_analysis(xa, weights, neofs=1, npcs=1, fn=None):
    """ Empirical Orthogonal Function of SST(t,x,y) field """
    assert type(xa)==xr.core.dataarray.DataArray
    assert type(weights)==xr.core.dataarray.DataArray
    assert 'time' in xa.dims
    assert np.shape(xa[0,:,:])==np.shape(weights)
    
    # anomalies by removing time mean
    xa = xa - xa.mean(dim='time')
    # Retrieve the leading EOF, expressed as the covariance between the leading PC
    # time series and the input xa anomalies at each grid point.
    solver = Eof(xa, weights=weights)
    eofs = solver.eofsAsCovariance(neofs=neofs)
    pcs = solver.pcs(npcs=npcs, pcscaling=1)
    if fn!=None:
        xr.merge([eofs, pcs]).to_netcdf(fn)
    return eofs, pcs


def IPO_TPI(run):
    """ raw timeseries of the three tripole SST indices described in Henley et al. (2015) """
    assert run in ['ctrl', 'rcp', 'lpd', 'lpi']
    if run in ['ctrl', 'rcp']:
        domain = 'ocn_rect'
    elif run in ['lpd', 'lpi']:
        domain = 'ocn_low'
    
    TPI1 = SST_index_from_monthly(run=run, index_loc=None, MASK=TPI_masks(domain, 1))
    TPI2 = SST_index_from_monthly(run=run, index_loc=None, MASK=TPI_masks(domain, 2))
    TPI3 = SST_index_from_monthly(run=run, index_loc=None, MASK=TPI_masks(domain, 3))
    
    TPI1.to_netcdf(f'{path_samoc}/SST/TPI1_monhtly_{run}.nc')
    TPI2.to_netcdf(f'{path_samoc}/SST/TPI2_monhtly_{run}.nc')
    TPI3.to_netcdf(f'{path_samoc}/SST/TPI3_monhtly_{run}.nc')
    
    return


def PMV_EOF_indices(run, extent):
    """ perform EOF of monthly, deseasonalized SST (detrended with global mean SST)
    NB: we will use 60S to 60N SST data as polar data is limited in observations

    1. for detrending: compute monthly global mean SST time series, deseasonalize them
        >>> see `SST_data_generation.py` file for scripts
    2. North Pacific monthly output fields
        2.1. create monthly SST field
             (if appropriate: determine extend of grid, limit all coordinates)
             save as single file
             a) North of 38 deg S
             b) North of Equator
             b) North of 20 deg N
        2.2. deseasonalize 
        2.3. detrend global mean, deseasonalized SST, 
        2.4. (remove mean at each point)
    3. EOF analysis
       --> index is first principal component
    4. regress time series on global maps
        goal: perform EOF of monthly, deseasonalized SST (detrended with global mean SST)
    NB: we will use 60S to 60N SST data as polar data is limited in observations

    
    
    4. regress time series on global maps
        >>> in correspoinding .ipynb files
    """
    assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'had']
    assert extent in ['38S', 'Eq', '20N']
    
    if run in ['ctrl', 'rcp']:
        domain = 'ocn_rect'
        run_name = f'rect_{run}'
    elif run in ['lpd', 'lpi']:
        domain = 'ocn_low'
        run_name = run
    elif run=='had':
        domain = 'ocn_had'
        run_name = run
    
    # load fields (generated in `SST_data_generatoin.py`)
    # <SST>(t):  60S-60N (tropical/extratropical) monthly timeseries
    SST_xm = xr.open_dataarray(f'{path_samoc}/SST/SST_60S_60N_mean_monthly_{run_name}.nc', decode_times=False)
    
    # SST(t,x,y): monthly SST field limited to Pacific North of `extent` (*)
    SST = xr.open_dataarray(f'{path_samoc}/SST/SST_monthly_{run_name}.nc', decode_times=False)
    
    print('opened datasets')
    
    # deseasonalize
    SST_xm_ds = deseasonalize(SST_xm )
    SST_ds    = deseasonalize(SST)
    
    print('deseasonalized datasets')
    
    # some of the time series are not the same length
    if run=='ctrl' :  SST_xm_ds = SST_xm_ds[:-7]
    elif run=='rcp':  SST_xm_ds = SST_xm_ds[:-1]
    
    # detrend
    SST_ds_dt = SST_ds - SST_xm_ds
    
    print('detrended SSTs')
    
    # remove mean at each point
    SST_ds_dt_dm = SST_ds_dt - SST_ds_dt.mean('time')
    
    # EOF analysis
    # N.B. cut off two year on either end as the arbitrary start month biases the filtered time series
    if extent=='38S':
        latS, lonE = -38, 300
    elif extent=='Eq':
        latS, lonE = 0, 285
    elif extent=='20N':
        latS, lonE = 20, 255
    
    AREA = xr_AREA(domain=domain)
    Pac_MASK = mask_box_in_region(domain=domain, mask_nr=2, bounding_lats=(latS,68), bounding_lons=(110,lonE))
    Pac_area = AREA.where(Pac_MASK)
    fn = f'{path_samoc}/SST/SST_PDO_EOF_{extent}_{run}.nc'
    
    print('prepared EOF')
    
    eof, pc  = EOF_SST_analysis(xa=SST_ds_dt_dm[24:-24,:,:].where(Pac_MASK),
                                weights=Pac_area, fn=fn)
    
    print('performed EOF')    
    
    return eof, pc


def bounding_lats_lons(index):
    """ bounding latitudes and longitudes """
    if index=='AMO':
        (blats, blons) = bll_AMO
        mask_nr = 6
    elif index=='SOM':
        (blats, blons) = bll_SOM
        mask_nr = 0
    elif index=='TPI1':  # use monthly data
        (blats, blons) = bll_TPI1
        mask_nr = 2
    elif index=='TPI2':
        (blats, blons) = bll_TPI3
        mask_nr = 2
    elif index=='TPI3':
        (blats, blons) = bll_TPI3
        mask_nr = 2
    return blats, blons, mask_nr
    


def SST_index(index, run, detrend_signal='GMST', time_slice='full'):
    """ calcalates SST time series from yearly detrended SST dataset """
    assert index in ['AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']
    assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'had']
    assert detrend_signal in ['GMST', 'AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']
    
    print(index, run)
    
    if run in ['ctrl', 'rcp']:
        domain = 'ocn'  #check this
        dims = ('nlat', 'nlon')
    elif run in ['lpd', 'lpi']:
        domain = 'ocn_low'
        dims = ('nlat', 'nlon')
    elif run=='had':
        domain = 'ocn_had'
        dims = ('latitude', 'longitude')
        
    
    blats, blons, mask_nr = bounding_lats_lons(index)
        
    MASK = mask_box_in_region(domain=domain, mask_nr=mask_nr, bounding_lats=blats, bounding_lons=blons)
    AREA = xr_AREA(domain=domain).where(MASK)
    index_area = AREA.sum()

    if run=='had' and detrend_signal in ['AMO', 'SOM'] or detrend_signal=='GMST':
        print(f'underlying SST field: detrended with {detrend_signal}, no filtering')
        if detrend_signal=='GMST': 
            print('GMST(t) signal scaled at each grid point\n')
        else:
            print(f'{detrend_signal}(t) removed from all SST gridpoints without scaling\n')
        if time_slice=='full':
            fn = f'{path_samoc}/SST/SST_{detrend_signal}_dt_yrly_{run}.nc'
            trange = ''
        else:
            first_year, last_year = determine_years_from_slice(run=run, tres='yrly', time_slice=time_slice)
            trange = f'_{first_year}_{last_year}'
            fn = f'{path_samoc}/SST/SST_{detrend_signal}_dt_yrly{trange}_{run}.nc'
        assert os.path.exists(fn)
        SST_yrly = xr.open_dataarray(fn).where(MASK)
        detr = f'_{detrend_signal}_dt'
        
    else:  # run=='had' and detrend_signal!='GMST'
        print('underlying SST field: no detrending, no filtering')
        if detrend_signal in ['AMO', 'SOM']:
            print(f'{detrend_signal} must subsequently be detrended with polynomial\n')
        else:
            print(f'{detrend_signal} must not be detrended since forcing signal compensated in TPI\n')
        SST_yrly = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_{run}.nc').where(MASK)
        detr = ''
        
    SSTindex = SST_area_average(xa_SST=SST_yrly, AREA=AREA, AREA_index=index_area, MASK=MASK, dims=dims)
    SSTindex.to_netcdf(f'{path_samoc}/SST/{index}{detr}_raw{trange}_{run}.nc')
    
    return SSTindex



def SST_remove_forced_signal(run, tres='yrly', detrend_signal='GMST', time_slice='full'):
    """ removed the scaled, forced GMST signal (method by Kajtar et al. (2019))
    
    1. load raw SST data
    2. generate forced signal (either quadtrend or CMIP MMEM)
    3. regress forced signal onto SST data -> \beta
    4. use regression coefficient \beta to generate SST signal due to forcing
    5. remove that signal
    
    run  ..
    tres .. time resolution
    detrend_signal .. either GMST (Kajtar et al. (2019))
                      or target region (Steinman et al. (2015))
    time_slice     .. time range selected
    """
    assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'had']
    assert tres in ['yrly', 'monthly']
    assert detrend_signal in ['GMST', 'AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']
    if detrend_signal in ['AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']:
        assert run=='had'
    if run=='had':
        assert time_slice=='full'
    
    


    # file name and domain
    fn = f'{path_samoc}/SST/SST_{tres}_{run}.nc'
    if run in ['ctrl', 'rcp']:
        if tres=='yrly':
            domain = 'ocn'
        elif tres=='monthly':
            domain = 'ocn_rect'
    elif run in ['lpd', 'lpi']:  
        domain = 'ocn_low'
    elif run=='had':
        domain = 'ocn_had'

    first_year, last_year = determine_years_from_slice(run, tres, time_slice)
    

#     sys.exit("Error message")
    # 1. load data
    MASK = boolean_mask(domain=domain, mask_nr=0, rounded=True)
    SST = xr.open_dataarray(f'{path_samoc}/SST/SST_{tres}_{run}.nc', decode_times=False).where(MASK)
    if time_slice is not 'full':
        assert type(time_slice)==tuple
        SST = SST.sel(time=slice(*time_slice))
        
    if tres=='monthly':  # deseasonalize
        for t in range(12):
            SST[t::12,:,:] -= SST[t::12,:,:].mean(dim='time')
        
    SST = SST - SST.mean(dim='time')
    
    # 2/3/4. calculate forced signal
    forced_signal = forcing_signal(run=run, tres=tres, detrend_signal=detrend_signal, time_slice=time_slice)
    
    if detrend_signal=='GMST':
        beta = lag_linregress_3D(forced_signal, SST)['slope']
        if run=='had':
            beta = xr.where(abs(beta)<5, beta, np.median(beta))
        ds = xr.merge([forced_signal, beta])
        ds.to_netcdf(f'{path_samoc}/SST/SST_beta_{detrend_signal}_{tres}_{run}.nc')
        forced_map = beta * forced_signal
        
        # 5.
        SST_dt = SST - forced_map
        SST_dt -= SST_dt.mean(dim='time')
        
    elif detrend_signal in ['AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']:
        # these indices will be detrended afterwards
        SST_dt = SST - forced_signal
        ds = None
        
    if time_slice=='full':
        SST_dt.to_netcdf(f'{path_samoc}/SST/SST_{detrend_signal}_dt_{tres}_{run}.nc')
    else:
        SST_dt.to_netcdf(f'{path_samoc}/SST/SST_{detrend_signal}_dt_{tres}_{first_year}_{last_year}_{run}.nc')
    
    return SST_dt, ds




def forcing_signal(run, tres, detrend_signal, time_slice='full'):
    """ GMST forced component
    run            .. dataset
    tres           .. time resolution
    detrend_signal .. 
    time_slice
    """
    assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'had']
    assert tres in ['yrly', 'monthly']
    assert detrend_signal in ['GMST', 'AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']
    
    if run in ['ctrl', 'rcp', 'lpd', 'lpi']:
        assert detrend_signal=='GMST'
        forced_signal = xr.open_dataset(f'{path_samoc}/GMST/GMST_{tres}_{run}.nc', decode_times=False).GMST
        if run=='rcp':
            forced_signal = xr_quadtrend(forced_signal)
        else:
            forced_signal = xr_lintrend(forced_signal)
        if tres=='yrly':
            times = forced_signal['time'] + 31 # time coordinates shifted by 31 days (SST saved end of January, GMST beginning)
            if run=='ctrl':  # for this run, sometimes 31 days, sometimes 15/16 days offset
                times = xr.open_dataset(f'{path_samoc}/SST/SST_yrly_ctrl.nc', decode_times=False).time
            forced_signal = forced_signal.assign_coords(time=times)
#         elif tres=='monthly':
            
    elif run=='had':
        forced_signal = xr.open_dataarray(f'{path_data}/CMIP5/KNMI_CMIP5_{detrend_signal}_{tres}.nc', decode_times=False)
        if tres=='monthly':  # deseasonalize
            for t in range(12):
                forced_signal[t::12] -= forced_signal[t::12].mean(dim='time')
        
        if tres=='yrly':# select 1870-2018
            times = (forced_signal['time'].astype(int) - 9)*365
            forced_signal = forced_signal.assign_coords(time=times)  # days since 1861
            forced_signal = forced_signal[9:158]
        elif tres=='monthly':
            # ...
            forced_signal = forced_signal[9*12:158*12-1]
            times = xr.open_dataarray(f'{path_samoc}/SST/SST_monthly_had.nc', decode_times=False).time.values
            forced_signal = forced_signal.assign_coords(time=times)
            
    if time_slice is not 'full':
        forced_signal = forced_signal.sel(time=slice(*time_slice))
    
    forced_signal -= forced_signal.mean()
    forced_signal.name = 'forcing'
    
    return forced_signal


def determine_years_from_slice(run, tres, time_slice):
    assert time_slice is not 'full'
    if tres=='yrly':
        first_year, last_year = int(time_slice[0]/365), int(time_slice[1]/365)
    elif tres=='monthly':
        if run in ['ctrl', 'rcp']:
            first_year, last_year = int(time_slice[0]/12), int(time_slice[1]/12)
            if run=='ctrl':
                first_year += 100
                last_year  += 100
            elif run=='rcp':
                first_year += 2000
                last_year  += 2000
        elif run in ['lpd', 'lpi']:
            first_year, last_year = int(time_slice[0]/365), int(time_slice[1]/365)
    return first_year, last_year



if __name__=="__main__":
    idx = sys.argv[1]
    run = sys.argv[2]
    
    assert idx in ['IPO']
    assert run in ['ctrl', 'rcp', 'lpd', 'lpi']
    
    print(f'running {idx} index of {run} run\n')
    print(f'started at\n{datetime.datetime.now()}\n\n')

    if idx=='IPO':
        IPO_TPI(run)
        
    print(f'\n\nfinished at\n{datetime.datetime.now()}')