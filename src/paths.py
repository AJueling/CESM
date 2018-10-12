import numpy as np

def CESM_filename(domain, run, y, m, name=None):
    """ filename creation 
    
    input:
    domain   .. (str) 'ocn' or 'atm'
    run .. (str) 'ctrl' or 'rcp'
    y        .. (int) year
    m        .. (int) month; if 0, then yearly file
    name     .. (str) added to yrly file
    
    output:
    file     .. (str) filename
    """
    assert type(y)==np.dtype(int) and type(m)==np.dtype(int)
    assert m<13
    
    time = f'{y:04}-{m:02}'
    
    if domain=='ocn':
        if run=='ctrl':
            if m==0:  # yearly files
                file = f'{path_yrly_ctrl}/ocn_yrly_{name}_{y:04}.nc'
            else:
                file = f'{path_ocn_ctrl}/{spinup}.pop.h.{time}.nc'
        elif run=='rcp':
            if m==0:
                file = f'{path_yrly_rcp}/ocn_yrly_{name}_{y:04}.nc'
            else:
                file = f'{path_ocn_rcp}/{rcpstr}.pop.h.{time}.nc'
                
    if domain=='atm':
        if run=='ctrl':
            if m==0:
                file = f'{path_yrly_ctrl}/atm_yrly_{name}_{y:04}.nc'
            else:
                file = f'{path_atm_ctrl}/{spinup}.cam2.h0.{time}.nc'
        elif run=='rcp':
            if m==0:
                file = f'{path_yrly_rcp}/atm_yrly_{name}_{y:04}.nc'
            else:
                file = f'{path_atm_rcp}/{rcpstr}.cam2.h0.{time}.nc'
            
    return file


# STRINGS

spinup = 'spinup_pd_maxcores_f05_t12'
rcpstr = 'rcp8.5_co2_f05_t12'


# PATHS

path_results = '/home/dijkbio/andre/CESM/results'

# CESM data
path_CESM = '/projects/0/prace_imau/prace_2013081679/cesm1_0_4'
path_samoc = '/projects/0/samoc/andre/CESM'

# grid
path_ocn_grid = f'{path_CESM}/inputdata/ocn/pop/tx0.1v2/grid/'

# currently running
path_run_ctrl = f'{path_CESM}/{spinup}/run'
path_run_rcp  = f'{path_CESM}/{rcpstr}/run'

# then copied to
path_ocn_ctrl = f'{path_CESM}/{spinup}/OUTPUT/ocn/hist/monthly'
path_atm_ctrl = f'{path_CESM}/{spinup}/OUTPUT/atm/hist/monthly'

path_ocn_rcp  = f'{path_CESM}/rcp8.5_co2_f05_t12/OUTPUT/ocn/hist/monthly'
path_atm_rcp  = f'{path_CESM}/rcp8.5_co2_f05_t12/OUTPUT/atm/hist/monthly'

# interpolated to rectangular 0.4 deg grid
path_ocn_ctrl_rect = f'{path_CESM}/{spinup}/OUTPUT/ocn/hist/monthly_rect'
path_atm_ctrl_rect = f'{path_CESM}/{spinup}/OUTPUT/atm/hist/monthly_rect'

path_ocn_rcp_rect  = f'{path_CESM}/{rcpstr}/OUTPUT/ocn/hist/monthly_rect'
path_atm_rcp_rect  = f'{path_CESM}/{rcpstr}/OUTPUT/atm/hist/monthly_rect'

# OHC files created by Rene
path_rene          = '/projects/0/prace_imau/prace_2013081679/rene/CESM'
path_ohc_rene      = f'{path_rene}/OceanHeatContent/Global_0.1'
path_ohc_rene_rect = f'{path_rene}/OceanHeatContent/Global_0.4'

# yearly files
path_yrly_ctrl = f'{path_samoc}/ctrl'
path_yrly_rcp  = f'{path_samoc}/rcp'


# FILES

# example files to use for tests
file_ex_ocn_ctrl = CESM_filename(domain='ocn', run='ctrl', y=200, m=1)
file_ex_ocn_rcp  = CESM_filename(domain='ocn', run='rcp', y=2000, m=1)
file_ex_atm_ctrl = CESM_filename(domain='atm', run='ctrl', y=200, m=1)

file_ex_ocn_rect  = f'{path_ocn_ctrl_rect}/{spinup}.pop.h.200-01.interp900x602.nc'

file_ex_ohc_hires = f'{path_ohc_rene}/OHC_0200-01_All.nc'

file_ex_ocn_TEMP_PD_yrly = f'{path_samoc}/ctrl/ocn_yrly_TEMP_PD_0200.nc'
file_ex_atm_T_T850_yrly  = f'{path_samoc}/ctrl/atm_yrly_T_T850_0200.nc'