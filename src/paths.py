# CESM data
path_CESM = '/projects/0/prace_imau/prace_2013081679/cesm1_0_4'

# grid
path_ocn_grid = f'{path_CESM}/inputdata/ocn/pop/tx0.1v2/grid/'

# currently running
path_run_ctrl = f'{path_CESM}/spinup_pd_maxcores_f05_t12/run'
path_run_rcp  = f'{path_CESM}/rcp8.5_co2_f05_t12/run'

# then copied to
path_ocn_ctrl = f'{path_CESM}/spinup_pd_maxcores_f05_t12/OUTPUT/ocn/hist/monthly'
path_ocn_rcp  = f'{path_CESM}/rcp8.5_co2_f05_t12/OUTPUT/ocn/hist/monthly'

# interpolated to rectangular 0.4 deg grid
path_ocn_ctrl_rect = f'{path_CESM}/spinup_pd_maxcores_f05_t12/OUTPUT/ocn/hist/monthly_rect'
path_ocn_rcp_rect  = f'{path_CESM}/rcp8.5_co2_f05_t12/OUTPUT/ocn/hist/monthly_rect'

# OHC files created by Rene
path_ohc_rene      = '/projects/0/prace_imau/prace_2013081679/rene/CESM/OceanHeatContent/Global_0.1'
path_ohc_rene_rect = '/projects/0/prace_imau/prace_2013081679/rene/CESM/OceanHeatContent/Global_0.4'