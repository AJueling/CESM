TO DO LIST
---

__16.01.2018__


_highest priority_
- OHC calculations for low res runs
- PDO timeseries to compare with TPI index


_middle priority_
- POP output to .nc
- OHC transport across basin boundaries
- AMWG (and other) package results
- MPI GE dataset

_lower priority_
- fix SHF calculations
- calculate SOM energy pathways
- OHC in 200-400 m in the SO for Sybren

_technical issues_
- code coverage (especially in xr_integrate)
- tests:
	- physical consistency
	- testing functions

__Literature__

- Mann paper regression of indices on OHC change

- apply same methodology to 



ideas:
- hiatuses in the future, both in CESM and CMIP5,
  maybe sea ice changes will change the OH uptake
- linear regression model of OHC changes on index timeseries
    1. index to dOHC (lagged correlation)
    2. combine to multivariate linear model
  maybe we can show that the skill improves when we incorporate the SOM
  (Mann et al. did a hiatus / indices model)
- regional change vs. global change (video of SST trends; safe havens)
- dOHC: forcing vs response under climate change

side projects:



integration functions:

xr_int_global(da, AREA, DZ)                            dx dy dz
xr_int_global_level(da, AREA, DZ)                      dx dy
xr_int_vertical(da, DZ)                                      dz
xr_int_zonal(da, HTN, LATS, AREA, DZ)                  dx    dz
xr_int_zonal_level(da, HTN, LATS, AREA, DZ, dx=1)      dx
xr_zonal_int_bins(da, LATS, AREA, dx=1)
lat_binning(dx)
xr_vol_int(xa, AREA, DZ, levels=False, zonal=False)
xr_int_along_axis(xa, DZ, axis)                        generic
xr_vol_int_regional(xa, AREA, DZ, MASK)
find_regional_coord_extent(MASK)
xr_vol_mean(xa, AREA, DZ)
xr_surf_int(xa, AREA)
xr_surf_mean(xa, AREA)
xr_zonal_mean(xa, AREA, dx, lat_name)


