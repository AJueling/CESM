# =============================================================================
# GENERATING SST DATA FILES
# =============================================================================

# %%time
# # ca. 4:30 min for ctrl/rcp, 1:25 for lpi
# # stacking files into one xr Dataset object
# for run in ['ctrl', 'rcp']:  # ['lpi', 'lpd']:
#     for i, (y,m,s) in enumerate(IterateOutputCESM('ocn', run, 'yrly', name='TEMP_PD')):
#         print(y)
#         da = xr.open_dataset(s, decode_times=False).TEMP[0,:,:]
#         da = da.drop(['z_t', 'ULONG', 'ULAT'])
#         da['TLAT' ] = da['TLAT' ].round(decimals=2)
#         da['TLONG'] = da['TLONG'].round(decimals=2)
#         del da.encoding["contiguous"]
#         ds = t2ds(da=da, name='SST', t=int(round(da.time.item())))
#         ds.to_netcdf(path=f'{path_samoc}/SST/SST_yrly_{run}_{y}.nc', mode='w')
     
#     combined = xr.open_mfdataset(f'{path_samoc}/SST/SST_yrly_{run}_*.nc',
#                                  concat_dim='time',
#                                  autoclose=True,
#                                  coords='minimal')
#     combined.to_netcdf(f'{path_samoc}/SST/SST_yrly_{run}.nc')
#     # remove extra netCDF files


# =============================================================================
# GLOBAL MEAN TIME SERIES
# =============================================================================

# %%time
# # 29 min
# TAREA = xr_AREA('ocn')
# MASK_ocn = boolean_mask('ocn', 0)
# global_area = TAREA.where(MASK_ocn).sum()
# for run in ['ctrl', 'rcp']:
#     for i, (y,m,s) in enumerate(IterateOutputCESM(domain='ocn', run=run, tavg='monthly')):
#         SST = xr.open_dataset(s, decode_times=False).TEMP[0,0,:,:]
#         SST_gm = SST_index(xa_SST=SST, AREA=TAREA, index_loc=global_ocean, AREA_index=global_area, MASK=MASK_ocn, dims=('nlat', 'nlon'))
#         if m==1: print(y, SST_gm.item()) 
#         if i==0:  SST_new = SST_gm
#         else:     SST_new = xr.concat([SST_new, SST_gm], dim='time')
#     SST_new.to_netcdf(f'{path_samoc}/SST/SST_global_mean_monthly_{run}.nc')

# %%time
# # CTRL and RCP global means: ocn_rect
# # 7 min for both
# AREA_rect = xr_AREA('ocn_rect')
# MASK_rect = boolean_mask('ocn_rect', 0)
# global_area2 = AREA_rect.where(MASK_rect).sum()
# for run in ['ctrl', 'rcp']:
#     for i, (y,m,s) in enumerate(IterateOutputCESM(domain='ocn_rect', run=run, tavg='monthly')):
#         SST = xr.open_dataset(s, decode_times=False).TEMP[0,:,:]
#         SST_gm = SST_index(xa_SST=SST, AREA=AREA_rect, index_loc=gl_ocean_rect,
#                            AREA_index=global_area2, MASK=MASK_rect,
#                            dims=('t_lat', 't_lon'))
#         if m==1: print(y, SST_gm.item()) 
#         if i==0:  SST_new = SST_gm
#         else:     SST_new = xr.concat([SST_new, SST_gm], dim='time')
#     SST_new.to_netcdf(f'{path_samoc}/SST/SST_global_mean_monthly_rect_{run}.nc')

# %%time
# # CTRL and RCP 60S to 60N means: ocn_rect
# # 7:26 min for both
# AREA_rect = xr_AREA('ocn_rect').sel({'t_lat':slice(-60,60)})
# MASK_rect = boolean_mask('ocn_rect', 0).sel({'t_lat':slice(-60,60)})
# global_area = AREA_rect.where(MASK_rect).sum()
# for run in ['ctrl', 'rcp']:
#     for i, (y,m,s) in enumerate(IterateOutputCESM(domain='ocn_rect', run=run, tavg='monthly')):
#         SST = xr.open_dataset(s, decode_times=False).TEMP[0,:,:].sel({'t_lat':slice(-60,60)})
#         SST_gm = SST_index(xa_SST=SST, AREA=AREA_rect, index_loc=gl_ocean_rect,
#                            AREA_index=global_area, MASK=MASK_rect,
#                            dims=('t_lat', 't_lon'))
#         if m==1: print(y, SST_gm.item()) 
#         if i==0:  SST_new = SST_gm
#         else:     SST_new = xr.concat([SST_new, SST_gm], dim='time')
#     SST_new.to_netcdf(f'{path_samoc}/SST/SST_60S_60N_mean_monthly_rect_{run}.nc')




# =============================================================================
# 60S-60N MEAN TIME SERIES
# =============================================================================

# %%time
# # LPI and LPD global means: ocn_low
# # 1 hr for LPD, 19 min for LPI
# AREA_low = xr_AREA('ocn_low')
# MASK_low = boolean_mask('ocn_low', 0)
# global_area_low = AREA_low.where(MASK_low).sum()
# for run in ['lpi']:#, 'lpd']:
#     for i, (y,m,s) in enumerate(IterateOutputCESM(domain='ocn_low', run=run, tavg='monthly')):
#         SST = xr.open_dataset(s, decode_times=False).TEMP[0,0,:,:]
#         SST_gm = SST_index(xa_SST=SST, AREA=AREA_low, index_loc=gl_ocean_low, 
#                            AREA_index=global_area_low, MASK=MASK_low, 
#                            dims=('nlat', 'nlon'))
#         if m==1: print(y, SST_gm.item()) 
#         if i==0:  SST_new = SST_gm
#         else:     SST_new = xr.concat([SST_new, SST_gm], dim='time')
#     SST_new.to_netcdf(f'{path_samoc}/SST/SST_global_mean_monthly_{run}.nc')



# HadISST 60S-60N time series is generated in SST_obs.ipynb

# =============================================================================
# ISOLATE PACIFIC SST DATA
# =============================================================================


# %%time
# # 8 mins for 200 years of ctrl
# n=1

# for j, run in enumerate(['ctrl', 'rcp']):
#     for y,m,s in IterateOutputCESM(domain='ocn_rect', tavg='monthly', run=run):
#         xa = xr.open_dataset(s, decode_times=False).TEMP[0,:,:]
#         if m==1:
#             print(y)
#             xa_out = xa.copy()    
#         else:
#             xa_out = xr.concat([xa_out, xa], dim='time')

#         if m==12:
#             xa_out.to_netcdf(f'{path_samoc}/SST/SST_monthly_rect_{run}_{y}.nc')
# run='ctrl'
# combined = xr.open_mfdataset(f'{path_samoc}/SST/SST_monthly_rect_{run}_*.nc', concat_dim='time', decode_times=False)
# combined.to_netcdf(f'{path_samoc}/SST/SST_monthly_rect_{run}.nc')
# combined.close()
# remove yearly files

# %%time
# # 6 mins for 200 years of ctrl
# n = 1
# for i, r in enumerate(['Pac_38S', 'Pac_Eq', 'Pac_20N']):
#     latS = [-38, 0, 20][i]
#     lonE = [300, 285, 255][i]
#     Pac_MASK = mask_box_in_region(domain='ocn_rect', mask_nr=2, bounding_lats=(latS,68), bounding_lons=(110,lonE))
#     Pac_MASK = Pac_MASK.where(Pac_MASK.t_lon+1/.6*Pac_MASK.t_lat<333,0)
#     NPac_area = xr_AREA('ocn_rect').where(Pac_MASK, drop=True)
#     for run  in ['ctrl', 'rcp']:
#         for y,m,s in IterateOutputCESM(domain='ocn_rect', tavg='monthly', run=run):
#             xa = xr.open_dataset(s, decode_times=False).TEMP[0,:,:].where(Pac_MASK, drop=True)
#             if m==1:
#                 print(y)
#                 xa_out = xa.copy()    
#             else:
#                 xa_out = xr.concat([xa_out, xa], dim='time')
#             if m==12:
#                 xa_out.to_netcdf(f'{path_samoc}/SST/SST_monthly_{r}_rect_{run}_{y}.nc')
#     combined = xr.open_mfdataset(f'{path_samoc}/SST/SST_monthly_{r}_rect_{run}_*.nc', concat_dim='time', decode_times=False)
#     combined.to_netcdf(f'{path_samoc}/SST/SST_monthly_{r}_rect_{run}.nc')
#     combined.close()



# %%time
# # 
# n = 1
# for i, r in enumerate(['Pac_38S', 'Pac_Eq', 'Pac_20N']):
#     latS = [-38, 0, 20][i]
#     lonE = [300, 285, 255][i]
#     Pac_MASK = mask_box_in_region(domain='ocn_low', mask_nr=2, bounding_lats=(latS,68), bounding_lons=(110,lonE))
#     NPac_area = xr_AREA('ocn_low').where(Pac_MASK, drop=True)
#     for run  in ['lpd', 'lpi']:
#         for y,m,s in IterateOutputCESM(domain='ocn_low', tavg='monthly', run=run):
#             xa = xr.open_dataset(s, decode_times=False).TEMP[0,:,:].where(Pac_MASK, drop=True)
#             if m==1:
#                 print(y)
#                 xa_out = xa.copy()    
#             else:
#                 xa_out = xr.concat([xa_out, xa], dim='time')
#             if m==12:
#                 xa_out.to_netcdf(f'{path_samoc}/SST/SST_monthly_{r}_{run}_{y}.nc')
#     combined = xr.open_mfdataset(f'{path_samoc}/SST/SST_monthly_{r}_{run}_*.nc', concat_dim='time', decode_times=False)
#     combined.to_netcdf(f'{path_samoc}/SST/SST_monthly_{r}_{run}.nc')
#     combined.close()
# # remove yearly files


# ds = xr.open_dataset(file_HadISST, decode_times=False)
# run = 'had'
# for i, r in enumerate(['Pac_38S', 'Pac_Eq', 'Pac_20N']):
#     latS = [-38, 0, 20][i]
#     lonE = [300, 285, 255][i]
#     Pac_MASK = mask_box_in_region(domain='ocn_had', mask_nr=2, bounding_lats=(latS,68), bounding_lons=(110,lonE))
#     da = ds.sst.where(Pac_MASK, drop=False)
#     da.to_netcdf(f'{path_samoc}/SST/SST_monthly_{r}_{run}.nc')






# =============================================================================
# remove yearly files
# =============================================================================

# def cleanup_yearly_files():
# return