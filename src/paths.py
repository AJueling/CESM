import numpy as np

# POP files: 
#  75-275: /projects/0/samoc/pop/tx0.1/output/run_henk_mixedbc/tavg_rectgrid/
# 276-326: /projects/0/samoc/pop/tx0.1/output/run_henk_mixedbc_extravars_viebahn/tavg_rectgrid/

def CESM_filename(domain, run, y, m, name=None):
    """ filename creation 
    
    input:
    domain   .. (str) 'ocn', 'ocn_rect', 'ocn_low', 'atm', 'ice'
    run      .. (str) 'ctrl', 'rcp', 'lpd', 'lpi', 'pop'
    y        .. (int) year
    m        .. (int) month; if 0, then yearly file
    name     .. (str) added to yrly file
    
    output:
    file     .. (str) filename
    """
    assert domain in ['ocn', 'ocn_rect', 'ocn_low', 'atm', 'ice']
    assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'pop']
    assert type(y)==np.dtype(int) and type(m)==np.dtype(int)
    assert m>=0 and m<13
    
    time = f'{y:04}-{m:02}'
    time2 = f'{y:04}{m:02}'
    
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
        elif run=='lpd':
            if m==0:
                file = f'{path_yrly_lpd}/ocn_yrly_{name}_{y:04}.nc'
            else:
                file = f'{path_ocn_lpd}/{lpdstr}.pop.h.{time}.nc'
        elif run=='lpi':
            if m==0:
                file = f'{path_yrly_lpi}/ocn_yrly_{name}_{y:04}.nc'
            else:
                file = f'{path_ocn_lpi}/{lpistr}.pop.h.{time}.nc'
                
    elif domain=='ocn_rect':
        if run=='ctrl':
            if m==0:  # yearly files
                file = f'{path_yrly_ctrl}/ocn_yrly_{name}_{y:04}.interp900x602.nc'
            else:
                file = f'{path_ocn_rect_ctrl}/{spinup}.pop.h.{time}.interp900x602.nc'
        elif run=='rcp':
            if m==0:
                file = f'{path_yrly_rcp}/ocn_yrly_{name}_{y:04}.interp900x602.nc'
            else:
                file = f'{path_ocn_rect_rcp}/{rcpstr}.pop.h.{time}.interp900x602.nc'
        elif run=='pop':
            if m==0:
                if y<276:
                    file = f'{path_ocn_rect_pop1}/yearly/{popstr}.avg{y:04}.nc'
                elif y>=276:
                    file = f'{path_ocn_rect_pop2}/{popstr}.avg{y:04}.nc'
            else:
                if y<276:
                    file = f'{path_ocn_rect_pop1}/monthly/{popstr}.{time2}.interp900x602.nc'
                elif y>=276:
                    file = f'{path_ocn_rect_pop2}/{popstr}.{time2}.interp900x602.nc'
                
    elif domain=='ocn_low':
        if run=='lpd':
            if m==0:  # yearly files
                file = f'{path_yrly_lpd}/ocn_yrly_{name}_{y:04}.nc'
            else:
                file = f'{path_ocn_lpd}/{lpdstr}.pop.h.{time}.nc'
        elif run=='lpi':
            if m==0:
                file = f'{path_yrly_lpi}/ocn_yrly_{name}_{y:04}.nc'
            else:
                file = f'{path_ocn_lpi}/{lpistr}.pop.h.{time}.nc'
    
    elif domain=='atm':
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
        elif run=='lpd':
            if m==0:
                file = f'{path_atm_lpd}/{lpdstr}.cam.h0.avg{y:04}.nc'
            else:
                raise ValueError('monthly files are not available for lpd run!')
        elif run=='lpi':
            if m==0:
                file = f'{path_yrly_lpi}/atm_yrly_{name}_{y:04}.nc'
            else:
                file = f'{path_atm_lpi}/{lpistr}.cam2.h0.{time}.nc'
                
    elif domain=='ice':
        if run=='ctrl':
            if m==0:
                file = f'{path_yrly_ctrl}/ice_yrly_{name}_{y:04}.nc'
            else:
                file = f'{path_ice_ctrl}/{spinup}.cice.h.{time}.nc'
        elif run=='rcp':
            if m==0:
                file = f'{path_yrly_rcp}/ice_yrly_{name}_{y:04}.nc'
            else:
                file = f'{path_ice_rcp}/{rcpstr}.cice.h.{time}.nc'
            
    return file


# STRINGS

spinup = 'spinup_pd_maxcores_f05_t12'
rcpstr = 'rcp8.5_co2_f05_t12'
lpdstr = 'spinup_B_2000_cam5_f09_g16'
lpistr = 'b.PI_1pic_f19g16_NESSC_control'
popstr = 't.t0.1_42l_nccs01'


# PATHS

# my output (write)
path_results = '/home/dijkbio/andre/CESM/results'
path_data    = '/home/dijkbio/andre/CESM/data'
path_samoc   = '/projects/0/samoc/andre/CESM'

# CESM data (should read only)
path_CESM = '/projects/0/prace_imau/prace_2013081679/cesm1_0_4'
path_ctrl = f'{path_CESM}/{spinup}'
path_rcp  = f'{path_CESM}/{rcpstr}'
path_lpd  = f'/projects/0/acc/cesm/cesm1_1_2/{lpdstr}'
path_lpi  = f'/projects/0/acc/cesm/cesm1_0_5/{lpistr}'
path_pop  = f'/projects/0/samoc/pop/tx0.1'

# grid
path_ocn_grid = f'{path_CESM}/inputdata/ocn/pop/tx0.1v2/grid/'

# currently running
path_run_ctrl = f'{path_ctrl}/run'
path_run_rcp  = f'{path_rcp}/run'

# then copied to
path_ocn_ctrl = f'{path_ctrl}/OUTPUT/ocn/hist/monthly'
path_atm_ctrl = f'{path_ctrl}/OUTPUT/atm/hist/monthly'
path_ice_ctrl = f'{path_ctrl}/OUTPUT/ice/hist/monthly'

path_ocn_rcp  = f'{path_rcp}/OUTPUT/ocn/hist/monthly'
path_atm_rcp  = f'{path_rcp}/OUTPUT/atm/hist/monthly'
path_ice_rcp  = f'{path_rcp}/OUTPUT/ice/hist/monthlies'

path_ocn_lpd  = f'{path_lpd}/OUTPUT/ocn/hist/monthly'
path_atm_lpd  = f'{path_lpd}/OUTPUT/atm/hist/yearly'
path_ice_lpd  = f'{path_lpd}/OUTPUT/ice/hist'

path_ocn_lpi  = f'{path_lpi}/OUTPUT/ocn/hist/monthly'
path_atm_lpi  = f'{path_lpi}/OUTPUT/atm/hist'


# interpolated to rectangular 0.4 deg grid
path_ocn_rect_ctrl = f'{path_ctrl}/OUTPUT/ocn/hist/monthly_rect'
path_atm_ctrl_rect = f'{path_ctrl}/OUTPUT/atm/hist/monthly_rect'

path_ocn_rect_rcp  = f'{path_rcp}/OUTPUT/ocn/hist/monthly_rect'
path_atm_rcp_rect  = f'{path_rcp}/OUTPUT/atm/hist/monthly_rect'

path_ocn_rect_pop1 = f'{path_pop}/output/run_henk_mixedbc/tavg_rectgrid'
path_ocn_rect_pop2 = f'{path_pop}/output/run_henk_mixedbc_extravars_viebahn/tavg_rectgrid'


# OHC files created by Rene
path_rene          = '/projects/0/prace_imau/prace_2013081679/rene/CESM'
path_ohc_rene      = f'{path_rene}/OceanHeatContent/Global_0.1'
path_ohc_rene_rect = f'{path_rene}/OceanHeatContent/Global_0.4'

# yearly files
path_yrly_ctrl = f'{path_samoc}/ctrl'
path_yrly_rcp  = f'{path_samoc}/rcp'
path_yrly_lpd  = f'{path_samoc}/lpd'
path_yrly_lpi  = f'{path_samoc}/lpi'

# FILES

grid_file  = f'{path_CESM}/inputdata/ocn/pop/tx0.1v2/grid/horiz_grid_200709.ieeer8'

# example files to use for tests
file_ex_ocn_ctrl = CESM_filename(domain='ocn', run='ctrl', y= 200, m=1)
file_ex_ocn_rcp  = CESM_filename(domain='ocn', run='rcp' , y=2000, m=1)
file_ex_ocn_lpd  = CESM_filename(domain='ocn', run='lpd' , y= 200, m=1)
file_ex_ocn_lpi  = CESM_filename(domain='ocn', run='lpi' , y=1600, m=1)

file_ex_ocn_rect  = f'{path_ocn_rect_ctrl}/{spinup}.pop.h.0200-01.interp900x602.nc'

file_ex_atm_ctrl = CESM_filename(domain='atm', run='ctrl', y= 200, m=1)
file_ex_atm_rcp  = CESM_filename(domain='atm', run='rcp' , y=2000, m=1)
file_ex_atm_lpd  = CESM_filename(domain='atm', run='lpd' , y= 200, m=0)
file_ex_atm_lpi  = CESM_filename(domain='atm', run='lpi' , y=3000, m=1)

file_ex_ice_rcp  = f'{path_ice_rcp}/{rcpstr}.cice.h.2000-01.nc'
file_ex_ice_yrly = f'{path_ice_rcp}/{rcpstr}.cice.h.avg2000.nc'

file_HadISST     = f'{path_data}/HadISST_sst.nc'

# derived data
file_ex_ohc_hires = f'{path_ohc_rene}/OHC_0200-01_All.nc'

file_ex_ocn_TEMP_PD_yrly = f'{path_samoc}/ctrl/ocn_yrly_TEMP_PD_0200.nc'
file_ex_atm_T_T850_U_V_yrly  = f'{path_samoc}/ctrl/atm_yrly_T_T850_U_V_0200.nc'


file_geometry = f'{path_ocn_grid}/dzbc_pbc_s2.0_200709.ieeer8'

file_RMASK_ocn      = f'{path_samoc}/grid/RMASK_ocn.nc'
file_RMASK_ocn_rect = f'{path_samoc}/grid/RMASK_ocn_rect.nc'
file_RMASK_ocn_low  = f'{path_samoc}/grid/RMASK_ocn_low.nc'
file_RMASK_ocn_had  = f'{path_samoc}/grid/RMASK_ocn_had.nc'