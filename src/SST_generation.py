""" ____ Regression Pattern Pipeline ____

input:
    run    .. simulations 'ctrl' & 'lpd', or observations 'had'

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

detrending methods:
    'had'          :  two-factor
    'ctrl' or 'lpd':  quadratic pointwise
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
idx = str(sys.argv[1])
assert run in ['ctrl', 'lpd', 'had']
assert idx in ['AMO', 'SOM', 'TPI', 'PMV']

def time_print():
    return f'\n   {datetime.datetime.now().time()}'

# ==================================================================================================
print(f'\nRegression Pattern Pipeline for {run} {idx}', time_print(), '\n')
# ==================================================================================================

# 10x 149 year time segments and one 250 year segment
times_ctrl = (list(zip(np.arange(51, 150, 10), np.arange(200, 301, 10))))
times_lpd  = (list(zip(np.arange(154, 253, 10), np.arange(303, 452, 10))))
times_ctrl.append((51,301))
times_lpd .append((154,404))

if run=='ctrl':  times = times_ctrl
elif run=='lpd':  times = times_lpd


# ==================================================================================================
print('1. monthly and yearly SST files', time_print())
# ==================================================================================================

# yearly data must be present for all index calculations
try:
    fn = f'{path_samoc}/SST/SST_yrly_{run}.nc'
    assert os.path.exists(fn)
    print(f'   file exists: {fn}')
    if run=='ctrl':
        fn = f'{path_prace}/SST/SST_yrly_rect_ctrl.nc'
        assert os.path.exists()
        print(f'   file exists: {fn}')
except:
    DS().generate_yrly_SST_files(run=run)
    
# monthly data necessary for TPI and PMV
if idx in ['TPI', 'PMV']:
    try:
        fn = f'{path_samoc}/SST/SST_monthly_{run}.nc'
        assert os.path.exists(fn)
        print(f'   file exists: {fn}')
        if run=='ctrl':
            fn = f'{path_prace}/SST/SST_yrly_rect_ctrl.nc'
            assert os.path.exists(fn)
            print(f'   file exists: {fn}')  # ocn_rect as well
    except:
        DS().generate_monthly_SST_files(run=run)
    

# ==================================================================================================
print('2. deseasonalize and detrended SST field', time_print())
# ==================================================================================================

# yrly data
if run=='had':  # two-factor detrending
    try:
        fn = f'{path_samoc}/SST/SST_GMST_tfdt_yrly_had.nc'
        assert os.path.exists(fn)
        print(f'   file exists: {fn}')
    except:
        DS().SST_remove_forced_signal(run='had', tres='yrly', detrend_signal='GMST', time_slice='full')

if run in ['ctrl', 'lpd']:  # quadratic detrending
    for time in times:
        try:
            fn = f'{path_samoc}/SST/SST_quadratic_pwdt_yrly_{run}_{time[0]}_{time[1]}.nc'
            assert os.path.exists(fn)
            print(f'   file exists: {fn}')
        except:
            DS().SST_pointwise_detrending(run=run, tres='yrly', degree=2, time_slice=tuple(time))

            
# monthly data
if idx in ['TPI', 'PMV']:
    if run=='had':
        
        try:  # deseasonalize
            assert os.path.exists(f'{path_prace}/SST/SST_monthly_deseasonalized_{runs}.nc')
        except:
            DS().deseasonalize_monthly_data(run='had')
            
        try:  # detrend
            assert os.path.exists(f'{path_prace}/SST/SST_monthly_ds_tfdt_had.nc')
        except:
            DS().detrend_monthly_obs_two_factor()
            
        if idx=='PMV':  # subselect Pacific data
            for extent in ['38S', 'Eq', '20N']:
                try:
                    fn = f'{path_samoc}/SST/SST_monthly_ds_dt_{extent}_had.nc'
                    assert os.path.exists(fn)
                except:
                    DS().isolate_Pacific_SSTs(run=run, extent=extent, time_slice=None)
        
        
    elif run in ['ctrl', 'lpd']:
        for time in times:          
            
            try:  # deseasonalize
                fn = f'{path_samoc}/SST/SST_monthly_deseasonalized_{run}_{time[0]}_{time[1]}.nc'
                assert os.path.exists(fn)
            except:
                DS().deseasonalize_monthly_data(run=run, time_slice=tuple(time))
                
            try:  # detrend
                fn = f'{path_samoc}/SST/SST_monthly_ds_dt_{run}_{time[0]}_{time[1]}.nc'
                assert os.path.exists(fn)
            except:
                DS().detrend_monthly_data_pointwise(run=run, time_slice=tuple(time))
    
            if idx=='PMV':  # subselect Pacific data
                for extent in ['38S', 'Eq', '20N']:
                    try:
                        fn = f'{path_samoc}/SST/SST_monthly_ds_dt_{extent}_{run}_{time[0]}_{time[1]}.nc'
                        assert os.path.exists(fn)
                    except:
                        DS().isolate_Pacific_SSTs(run=run, extent=extent, time_slice=tuple(time))
        
        
# ==================================================================================================
print('3. raw SST indices', time_print())
# ==================================================================================================

# yrly data
if idx in ['AMO', 'SOM']:    dt = 'dt'
elif idx in ['TPI', 'PMV']:  dt = 'ds_dt'

if idx in ['AMO', 'SOM', 'TPI']:  # area average SST indices
    if run=='had':
        try:
            fn = f'{path_samoc}/SST/{idx}_{dt}_raw_{run}.nc'
            assert os.path.exists(fn)
        except:
            AI().derive_SST_avg_index(run=run, index=idx, time=None)

    elif run in ['ctrl', 'lpd']:
        for time in times:
            try:
                fn = f'{path_samoc}/SST/{idx}_{dt}_raw_{run}_{time[0]}_{time[1]}.nc'
                assert os.path.exists(fn)
            except:
                AI().derive_SST_avg_index(run=run, index=idx, time=time)
        
if idx=='PMV':  # EOF analysis
    for extent in ['38S', 'Eq', '20N']:
        if run=='had':
            try:
                fn = f'{path_samoc}/SST/PMV_EOF_{extent}_{run}.nc'
                assert os.path.exists(fn)
            except:
                AI().Pacific_EOF_analysis(run=run, extent=extent, time=None)
        elif run in ['ctrl', 'lpd']:
            for time in times:
                try:
                    fn = f'{path_samoc}/SST/PMV_EOF_{extent}_{run}_{time[0]}_{time[1]}.nc'
                    assert os.path.exists(fn)
                except:
                    AI().Pacific_EOF_analysis(run=run, extent=extent, time=time)

                    
# ==================================================================================================
print('4. final indices', time_print())
# ==================================================================================================

# make final indices
if run=='had':
    try:
        if idx=='PMV':
            for extent in ['38S', 'Eq', '20N']:
                fn = f'{path_samoc}/SST/PMV_EOF_{extent}_{run}.nc'
                assert os.path.exists(fn)
        else:
            fn = f'{path_samoc}/SST/{idx}_{run}.nc'
            assert os.path.exists(fn)
    except:
        AI().derive_final_SST_indices(self, run, index, time=None)

elif run in ['ctrl', 'lpd']:
    for time in times:
        try:
            if idx=='PMV':
                for extent in ['38S', 'Eq', '20N']:
                    fn = f'{path_samoc}/SST/PMV_EOF_{extent}_{run}.nc'
                    assert os.path.exists(fn)
            else:
                fn = f'{path_samoc}/SST/{idx}_{run}_{time[0]}_{time[1]}.nc'
                assert os.path.exists(fn)
        except:
            AI().derive_final_SST_indices(self, run, index, time=time)

        
# ==================================================================================================
print('5. autocorrelation fields', time_print())
# ==================================================================================================

for tavg in ['yrly', 'monthly']:
    if run=='had':
        try:
            fn = f'{path_samoc}/SST/SST_{tavg}_autocorrelation_{run}.nc'
            assert os.path.exists(fn)
        except:
            AI().derive_yrly_autocorrelations(run=run, tavg=tavg, time=None)

    elif run in ['ctrl', 'lpd']:
        for time in times:
            try:
                fn = f'{path_samoc}/SST/SST_{tavg}_autocorrelation_{run}_{time[0]}_{time[1]}.nc'
                assert os.path.exists(fn)
            except:
                AI().derive_yrly_autocorrelations(run=run, tavg=tavg, time=time)


# ==================================================================================================
print('6. SST regression on indices', time_print())
# ==================================================================================================

if run=='had':
    try:
        if idx=='PMV':  # EOF analysis
            for extent in ['38S', 'Eq', '20N']:
        else:
            fn = f'{path_samoc}/SST/{idx}_regr_{run}.nc'
            assert os.path.exists(fn)
    except:
        AI().make_yrly_regression_files(run, time=None)
        
elif run in ['ctrl', 'lpd']:
    for time in times:
        try:
            if idx=='PMV':  # EOF analysis
                for extent in ['38S', 'Eq', '20N']:
                    fn = f''
            else:
                fn = f'{path_samoc}/SST/{idx}_regr_{run}.nc'
                assert os.path.exists(fn)
        except:
            AI().make_yrly_regression_files(run, time=None)
        

        

# ==================================================================================================
print('6. time segments', time_print())
# ==================================================================================================
        
# time segments
if run in ['ctrl', 'lpd']:
    for t in starts:
        tslice = (t, t+148)
        
        print('   deriving final indices')
        try:
            for idx in ['AMO', 'SOM', 'TPI']:
                fn = f'{path_samoc}/SST/{idx}_quadratic_pwdt_{run}_{tslice[0]}_{tslice[1]}.nc'
                assert os.path.exists(fn)
            print(f'   filtered index files for {run} of segment {tslice} exist')
        except:
            AI().derive_final_SST_indices(run=run, dts='quadratic', tslice=tslice)
        
        print('   deriving autocorrelations')
        try:
            fn = f'{path_samoc}/SST/SST_autocorrelation_{run}_{tslice[0]}_{tslice[1]}.nc'
            assert os.path.exists(fn)
            print(f'   file exists: {fn}')
        except:
            AI().derive_yrly_autocorrelations(run, tslice)
        
        print('   deriving regression fields')
        try:
            for idx in ['AMO', 'SOM', 'TPI']:
                fn = f'{path_samoc}/SST/{idx}_regr_{run}_{tslice[0]}_{tslice[1]}.nc'
                assert os.path.exists(fn)
            print(f'   regression files for {run} of segment {tslice} exist')
        except:
            AI().make_yrly_regression_files(run, tslice)


# ==================================================================================================
print('\nSuccess.', time_print(), '\n')
# ==================================================================================================