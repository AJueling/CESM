""" ____ Regression Pattern Pipeline ____

input:
    run    .. simulations 'ctrl' & 'lpd', or observations 'had'
    idx    .. index 'AMO', 'SOM', 'TPI', 'PMV'

procedure:
    1. creating monthly and yearly SST files
    2. deseasonalizing & detrending of SST field
        - model data: pointwise quadratic fit
        - observations: two-factor detrending
    3. derive raw SST indices
    4. derive final SST indices

detrending methods:
    'had'          :  two-factor => anthro & natural CMIP5 MMM
    'ctrl' or 'lpd':  quadratic pointwise
"""

import os
import sys
import datetime
import numpy as np

from paths import path_prace
from ab_derivation_SST import DeriveSST as DS
from bc_analysis_fields import AnalyzeField as AF
from bd_analysis_indices import AnalyzeIndex as AI

def time_print():  return f'\n   {datetime.datetime.now().time()}'

def trex(fn, fct, kwargs={}):
    """  try: existence / except: create
    checks whether file exists already and if not, creates it
    input:
    fn      .. file name of file that function is supposed to create
    fct     .. function that will create data stored as `fn`
    kwargs  .. keyword arguments for function `fct`
    """
    assert type(fn)==str 
    assert hasattr(fct, '__call__')
    assert type(kwargs)==dict
    try:
        assert os.path.exists(fn), f'   now creating file:  {fn[46:]}'
        print(f'   file exists:  {fn[46:]}')
    except:
        fct(**kwargs)
    return

# 10x 149 year time segments and one 250 year segment
times_ctrl = list(zip(np.append(np.arange( 51, 150, 10),[ 51]),
                      np.append(np.arange(200, 299, 10),[301])))
times_lpd  = list(zip(np.append(np.arange(154, 253, 10),[154]),
                      np.append(np.arange(303, 402, 10),[404])))
times_had  = [None]

if __name__=='__main__':
    run = str(sys.argv[1])
    idx = str(sys.argv[2])
    assert run in ['ctrl', 'lpd', 'had']
    assert idx in ['AMO', 'SOM', 'TPI', 'PMV']

    # ==============================================================================================
    print(f'\nRegression Pattern Pipeline for {run} {idx}', time_print(), '\n')
    # ==============================================================================================

    if run=='ctrl':  times = times_ctrl
    elif run=='lpd':  times = times_lpd
    elif run=='had':  times = times_had

    # ==============================================================================================
    print('1. create monthly and yearly SST files', time_print())
    # ==============================================================================================

    kwargs = dict(run=run)

    # yearly data must be present for all index calculations
    fn = f'{path_prace}/SST/SST_yrly_{run}.nc'
    fct = DS().generate_yrly_SST_files
    trex(fn=fn, fct=fct, kwargs=kwargs)

    # ocn rect 
    if run=='ctrl':
        fn = f'{path_prace}/SST/SST_yrly_rect_ctrl.nc'
        trex(fn=fn, fct=fct, kwargs=kwargs)
        
    # monthly data necessary for TPI and PMV
    if idx in ['TPI', 'PMV']:
        fn = f'{path_prace}/SST/SST_monthly_{run}.nc'
        fct = DS().generate_monthly_SST_files
        trex(fn=fn, fct=fct, kwargs=kwargs)


    # ==============================================================================================
    print('\n2. deseasonalize and detrended SST field', time_print())
    # ==============================================================================================

    # yrly data
    for time in times:
        if run=='had':  # two-factor detrending
            fn = f'{path_prace}/SST/SST_yrly_tfdt_had.nc'
            fct = DS().SST_remove_forced_signal
            kwargs = dict(run='had', tavg='yrly', detrend_signal='GMST', time=None)
        elif run in ['ctrl', 'lpd']:  # quadratic detrending
            fn = f'{path_prace}/SST/SST_yrly_pwdt_{run}_{time[0]}_{time[1]}.nc'
            fct = DS().SST_pointwise_detrending
            kwargs = dict(run=run, tavg='yrly', degree=2, time=time)
        trex(fn=fn, fct=fct, kwargs=kwargs)

    # monthly data
    if idx in ['TPI', 'PMV']:
        for time in times:
            # deseasonalize
            fct = DS().deseasonalize_monthly_data
            if run=='had':
                fn = f'{path_prace}/SST/SST_monthly_ds_had.nc'
                kwargs = dict(run='had')
            elif run in ['ctrl', 'lpd']:
                fn = f'{path_prace}/SST/SST_monthly_ds_{run}_{time[0]}_{time[1]}.nc'
                kwargs = dict(run=run, time=tuple(time))
            trex(fn=fn, fct=fct, kwargs=kwargs)
                    
            # detrend
            if run=='had':
                fn = f'{path_prace}/SST/SST_monthly_ds_dt_had.nc'
                fct = DS().detrend_monthly_obs_two_factor
                kwargs = {}
            elif run in ['ctrl', 'lpd']:
                fn = f'{path_prace}/SST/SST_monthly_ds_dt_{run}_{time[0]}_{time[1]}.nc'
                fct = DS().detrend_monthly_data_pointwise
                kwargs = dict(run=run, time=tuple(time))
            trex(fn=fn, fct=fct, kwargs=kwargs)
            
            # subselect Pacific data
            if idx=='PMV':
                fct = DS().isolate_Pacific_SSTs
                for extent in ['38S', 'Eq', '20N']:
                    kwargs = dict(run=run, extent=extent, time=time)
                    if run=='had':
                        fn = f'{path_prace}/SST/SST_monthly_ds_dt_{extent}_had.nc'
                    elif run in ['ctrl', 'lpd']:
                        fn = f'{path_prace}/SST/SST_monthly_ds_dt_{extent}_{run}_{time[0]}_{time[1]}.nc'
                    trex(fn=fn, fct=fct, kwargs=kwargs)    

            
    # ==============================================================================================
    print('\n3. raw SST indices', time_print())
    # ==============================================================================================

    if idx in ['AMO', 'SOM']:    dt = 'dt'
    elif idx in ['TPI', 'PMV']:  dt = 'ds_dt'

    # area average SST indices
    if idx in ['AMO', 'SOM', 'TPI']:
        fct = AI().derive_SST_avg_index
        for time in times:
            kwargs = dict(run=run, index=idx, time=time)
            if run=='had':
                fn = f'{path_prace}/SST/{idx}_{dt}_raw_had.nc'
            elif run in ['ctrl', 'lpd']:
                fn = f'{path_prace}/SST/{idx}_{dt}_raw_{run}_{time[0]}_{time[1]}.nc'
            trex(fn=fn, fct=fct, kwargs=kwargs) 

    # EOF analysis
    if idx=='PMV':
        fct = AI().Pacific_EOF_analysis
        for extent in ['38S', 'Eq', '20N']:
            for time in times:
                kwargs = dict(run=run, extent=extent, time=time)
                if run=='had':
                    fn = f'{path_prace}/SST/PMV_EOF_{extent}_had.nc'
                elif run in ['ctrl', 'lpd']:
                    fn = f'{path_prace}/SST/PMV_EOF_{extent}_{run}_{time[0]}_{time[1]}.nc'
                trex(fn=fn, fct=fct, kwargs=kwargs)

                     
    # ==============================================================================================
    print('\n4. final indices', time_print())
    # ==============================================================================================

    fct = AI().derive_final_SST_indices

    for time in times:
        kwargs = dict(run=run, index=idx, time=time)
        if idx=='PMV':
            for extent in ['38S', 'Eq', '20N']:
                if run=='had':
                    fn = f'{path_prace}/SST/PMV_EOF_{extent}_{run}.nc'
                elif run in ['ctrl', 'lpd']:
                    fn = f'{path_prace}/SST/PMV_EOF_{extent}_{run}_{time[0]}_{time[1]}.nc'
        else:
            if run=='had':
                fn = f'{path_prace}/SST/{idx}_{run}.nc'
            elif run in ['ctrl', 'lpd']:
                fn = f'{path_prace}/SST/{idx}_{run}_{time[0]}_{time[1]}.nc'
        trex(fn=fn, fct=fct, kwargs=kwargs)


    # ==============================================================================================
    print('\n5. autocorrelation fields', time_print())
    # ==============================================================================================

    fct = AI().derive_autocorrelation_maps

    if idx in ['AMO', 'SOM']:  tavg = 'yrly'
    if idx in ['TPI', 'PMV']:  tavg = 'monthly'

    for time in times:
        kwargs = dict(run=run, tavg=tavg, time=time)
        if run=='had':
            fn = f'{path_prace}/SST/SST_{tavg}_autocorrelation_had.nc'
        elif run in ['ctrl', 'lpd']:
            fn = f'{path_prace}/SST/SST_{tavg}_autocorrelation_{run}_{time[0]}_{time[1]}.nc'
        trex(fn=fn, fct=fct, kwargs=kwargs)



    # ==============================================================================================
    print('6. SST regression on indices', time_print())
    # ==============================================================================================

    # fct = AI().make_regression_files

    # for time in times:
    #     kwargs = dict(run=run, index=idx, time=None)
    #     if idx=='PMV':
    #         for extent in ['38S', 'Eq', '20N']:
    #             if run=='had':
    #                 fn = f'{path_prace}/SST/{idx}_{extent}_regr_{run}.nc'
    #             elif run in ['ctrl', 'lpd']:
    #                 fn = f'{path_prace}/SST/{idx}_{extent}_regr_{run}_{time[0]}_{time[1]}.nc'
    #     else:
    #         if run=='had':
    #             fn = f'{path_prace}/SST/{idx}_regr_{run}.nc'
    #         elif run in ['ctrl', 'lpd']:
    #             fn = f'{path_prace}/SST/{idx}_regr_{run}_{time[0]}_{time[1]}.nc'
    #     trex(fn=fn, fct=fct, kwargs=kwargs)



    # ==============================================================================================
    print('7. spectra', time_print())
    # ==============================================================================================




    # ==============================================================================================
    print('\nSuccess.', time_print(), '\n')
    # ==============================================================================================