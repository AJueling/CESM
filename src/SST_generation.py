""" ____ Regression Pattern Pipeline ____

input:
    run    .. simulations 'ctrl' & 'lpd', or observations 'had'
    tslice .. timeslice: either all years 'full' or segment '(start year, end year)'

procedure:
    1. creating yearly SST files
    2. detrending of SST field
        a. segment-wise linear
        b. pointwise fit to GMST/GMSST quadratic fit
        c. pointwise quadratic fit
    3. derive raw SST indices
    4. derive final SST indices
    5. autocorrelation fields
    6. regression of SST fields on indices
    7. regression patterns: plots, quatify difference to observations
    8. stationarity: test pattern stbility through time

detrending methods:
    'had': two-factor
    'ctrl' or 'lpd': quadratic pointwise
"""

import os
import sys
import datetime
import numpy as np

from paths import path_samoc
from ab_derivation_SST import DeriveSST as DS
from bc_analysis_fields import AnalyzeField as AF
from bd_analysis_indices import AnalyzeIndex as AI

run = str(sys.argv[1])

assert run in ['ctrl', 'lpd', 'had']

def time_print():
    return f'\n   {datetime.datetime.now().time()}'


# ==================================================================================================
print(f'Regression Pattern Pipeline for {run}', time_print(), '\n')
# ==================================================================================================

if run=='ctrl':   starts = np.arange(50, 151, 5)
elif run=='lpd':  starts = np.arange(154, 415, 10)


# ==================================================================================================
print('1. yearly SST files', time_print())
# ==================================================================================================

try:
    fn = f'{path_samoc}/SST/SST_yrly_{run}.nc'
    assert os.path.exists(fn)
    print(f'   file exists: {fn}')
except:
    if run in ['ctrl', 'lpd']:  # for 'had' this is done elsewhere
        DS().generate_yrly_SST_files(run)


# ==================================================================================================
print('2. detrended SST field', time_print())
# ==================================================================================================
#     a. segment-wise linear
#     b. pointwise fit to GMST/GMSST quadratic fit
#     c. pointwise quadratic fit

if run=='had':
    try:
        for dt in ['sfdt', 'tfdt']:
            fn = f'{path_samoc}/SST/SST_GMST_{dt}_yrly_had.nc'
            assert os.path.exists(fn)
            print(f'   file exists: {fn}')
    except:
        DS().SST_remove_forced_signal(run='had', tres='yrly', detrend_signal='GMST', time_slice='full')

if run in ['ctrl', 'lpd']:
    if run=='ctrl':
        tslice = (50, 301)
        print(f'   for `ctrl` run only use years {tslice} for quadratic detrending')
    elif run=='lpd':
        tslice = 'full'

    # pointwise quadratic detrending
    try:
        fn = f'{path_samoc}/SST/SST_quadratic_pwdt_yrly_{run}.nc'
        assert os.path.exists(fn)
        print(f'   file exists: {fn}')
    except:
        DS().SST_pointwise_detrending(run=run, tres='yrly', degree=2, time_slice=tslice)

    # remove scaled GMST/GMSST fit, does not work ...

#     try:
#         fn = f'{path_samoc}/SST/GMSST_yrly_{run}.nc'
#         assert os.path.exists(fn)
#         print(f'file exists: {fn}')
#     except:
#         DS().generate_yrly_global_mean_SST(run=run)
#
#     try:
#         fn = f'{path_samoc}/SST/SST_GMST_sqdt_yrly_{run}.nc'
#         assert os.path.exists(fn)
#         print(f'  file exists: {fn}')
#     except:
#         DS().SST_remove_forced_signal(run=run, tres='yrly', detrend_signal='GMST', time_slice='full')


# ==================================================================================================
print('3. raw SST indices', time_print())
# ==================================================================================================

if run=='had':
    try:
        for idx in ['AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']:
            fn = f'{path_samoc}/SST/{idx}_GMST_tfdt_raw_{run}.nc'
            assert os.path.exists(fn)
        print(f'   raw index files for {run} exist')
    except:
        AI().derive_all_SST_avg_indices(run, dts='GMST', tslice='full')

if run in ['ctrl', 'lpd']:


    # quadratic pointwise
    try:
        for idx in ['AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']:
            fn = f'{path_samoc}/SST/{idx}_quadratic_pwdt_raw_{run}.nc'
            assert os.path.exists(fn)
        print(f'   raw index files for {run} exist')
    except:
        AI().derive_all_SST_avg_indices(run, dts='quadratic', tslice='full')

#     # quadratic GMST/GMSST fit
#     try:
#         for idx in ['AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']:
#             fn = f'{path_samoc}/SST/{idx}_GMST_sqdt_raw_{run}.nc'
#             assert os.path.exists(fn)
#         print(f'raw index files for {run} exist')
#     except:
#         AI().derive_all_SST_avg_indices(run, dts='GMST', tslice='full')


# ==================================================================================================
print('4. final SST indices', time_print())
# ==================================================================================================

if run=='had':
    AI().derive_final_SST_indices(run=run, dts='GMST', tslice='full')

if run in ['ctrl', 'lpd']:

    # full time series
    try:
        for idx in ['AMO', 'SOM', 'TPI']:
            fn = f'{path_samoc}/SST/{idx}_quadratic_pwdt_{run}_{tslice[0]}_{tslice[1]}.nc'
            assert os.path.exists(fn)
        print(f'   filtered index files for {run} for full time series exist')
    except:
        AI().derive_final_SST_indices(run=run, dts='quadratic', tslice='full')
#         AI().derive_final_SST_indices(run=run, dts='GMST', tslice='full')

    # time segments
    for t in starts:
        tslice = (t, t+148)
        try:
            for idx in ['AMO', 'SOM', 'TPI']:
                fn = f'{path_samoc}/SST/{idx}_quadratic_pwdt_{run}_{tslice[0]}_{tslice[1]}.nc'
                assert os.path.exists(fn)
            print(f'   filtered index files for {run} of segment {tslice} exist')
        except:
            AI().derive_final_SST_indices(run=run, dts='quadratic', tslice=tslice)


# ==================================================================================================
print('5. autocorrelation fields', time_print())
# ==================================================================================================

# for full fields
try:
    fn = f'{path_samoc}/SST/SST_autocorrelation_{run}.nc'
    assert os.path.exists(fn)
    print(f'   file exists: {fn}')
except:
    AI().derive_yrly_autocorrelations(run, 'full')

# time segments
if run in ['ctrl', 'lpd']:
    for t in starts:
        tslice = (t, t+148)
        try:
            fn = f'{path_samoc}/SST/SST_autocorrelation_{run}_{tslice[0]}_{tslice[1]}.nc'
            assert os.path.exists(fn)
            print(f'   file exists: {fn}')
        except:
            AI().derive_yrly_autocorrelations(run, tslice)


# ==================================================================================================
print('6. SST regression on indices', time_print())
# ==================================================================================================

if run=='had':
    try:
        for idx in ['AMO', 'SOM', 'TPI']:
            fn = f'{path_samoc}/SST/{idx}_regr_{run}.nc'
            assert os.path.exists(fn)
        print(f'   regression files for {run} exist')
    except:
        AI().make_yrly_regression_files(run, 'full')

# time segments
if run in ['ctrl', 'lpd']:
    for t in starts:
        tslice = (t, t+148)
        try:
            for idx in ['AMO', 'SOM', 'TPI']:
                fn = f'{path_samoc}/SST/{idx}_regr_{run}_{tslice[0]}_{tslice[1]}.nc'
                assert os.path.exists(fn)
            print(f'   regression files for {run} of segment {tslice} exist')
        except:
            AI().make_yrly_regression_files(run, tslice)


# ==================================================================================================
print('\nSuccess.', time_print())
# ==================================================================================================