# CESM data
path_CESM = '/projects/0/prace_imau/prace_2013081679/cesm1_0_4'

# grid
path_ocn_grid = f'{path_CESM}/inputdata/ocn/pop/tx0.1v2/grid/'

# currently running
path_run_ctrl = f'{path_CESM}/spinup_pd_maxcores_f05_t12/run'
path_run_rcp  = f'{path_CESM}/rcp8.5_co2_f05_t12/run'

# then copied to
path_ocn_ctrl = f'{path_CESM}/spinup_pd_maxcores_f05_t12/OUTPUT/ocn/hist/monthly'
path_atm_ctrl = f'{path_CESM}/spinup_pd_maxcores_f05_t12/OUTPUT/atm/hist/monthly'

path_ocn_rcp  = f'{path_CESM}/rcp8.5_co2_f05_t12/OUTPUT/ocn/hist/monthly'
path_atm_rcp  = f'{path_CESM}/rcp8.5_co2_f05_t12/OUTPUT/atm/hist/monthly'

# interpolated to rectangular 0.4 deg grid
path_ocn_ctrl_rect = f'{path_CESM}/spinup_pd_maxcores_f05_t12/OUTPUT/ocn/hist/monthly_rect'
path_atm_ctrl_rect = f'{path_CESM}/spinup_pd_maxcores_f05_t12/OUTPUT/atm/hist/monthly_rect'

path_ocn_rcp_rect  = f'{path_CESM}/rcp8.5_co2_f05_t12/OUTPUT/ocn/hist/monthly_rect'
path_atm_rcp_rect  = f'{path_CESM}/rcp8.5_co2_f05_t12/OUTPUT/atm/hist/monthly_rect'

# OHC files created by Rene
path_ohc_rene      = '/projects/0/prace_imau/prace_2013081679/rene/CESM/OceanHeatContent/Global_0.1'
path_ohc_rene_rect = '/projects/0/prace_imau/prace_2013081679/rene/CESM/OceanHeatContent/Global_0.4'

# example files to use for tests
file_ex_ocn_hires = f'{path_ocn_ctrl}/spinup_pd_maxcores_f05_t12.pop.h.0200-01.nc' 
file_ex_atm_hires = f'{path_atm_ctrl}/spinup_pd_maxcores_f05_t12.cam2.h0.0200-01.nc'

file_ex_ocn_rect  = f'{path_ocn_ctrl_rect}/spinup_pd_maxcores_f05_t12.pop.h.200-01.interp900x602.nc'

file_ex_ohc_hires = f'{path_ohc_rene}/OHC_0200-01_All.nc'