""" Meridional Freshwater/Salt Transport

calculate the freshwater and salt transport terms
(overturning, azonal, total, and eddy)

"""
import os
import sys
import tqdm
import xgcm
import dask
import numpy as np
import xarray as xr
import pop_tools

from paths import file_ex_ocn_ctrl, file_ex_ocn_lpd, path_prace
from constants import rho_sw
from datetime import datetime
from xr_DataArrays import xr_DZ_xgcm

metrics = {
    ('X'): ['DXT', 'DXU', 'HTN', 'HUS'],  # X distances
    ('Y'): ['DYT', 'DYU', 'HTE', 'HUW'],  # Y distances
    ('Z'): ['DZT', 'DZU'],                # Z distances
    ('X', 'Y'): ['TAREA', 'UAREA']        # Areas
}
coords = {
    'X': {'center':'nlon_t', 'right':'nlon_u'},
    'Y': {'center':'nlat_t', 'right':'nlat_u'},
    'Z': {'center':'z_t', 'left':'z_w_top', 'right':'z_w_bot'}
}


def calc_transports_pbc(grid, ds, MASK_tt, MASK_uu, S0=35.0):
    """ calculate the merid. FW and SALT transport terms with xgcm grid including PBCs """
    # SALT
    SALT = grid.interp(grid.interp(ds.SALT.where(MASK_tt).where(ds.SALT>0), 'Y'), 'X')  # SALT interpolated to UU grid
    VVEL = (ds.VVEL*ds.DZU/ds.dz).where(MASK_uu).where(ds.VVEL<1000)  # reduced velocity due to PBCs

    # non-PBC-weighted mean here, as salinity differences between neighbouring cells are assumed very small
    SALT_xmean = SALT.mean('nlon_u')
    SALT_xmean.name = 'SALT_xmean'
    
    # VVEL terms:  weighted by cell thickness
    VVEL_xint  = (VVEL*ds.DXU).sum('nlon_u')  # zonal integral velocity (y,z) [cm^2/s]
    VVEL_xmean = (VVEL*ds.DXU).sum('nlon_u')/ds.DXU.where(VVEL<1000).sum('nlon_u')  # zonal mean velocity     (y,z) [cm/s]
    VVEL_xmean.name = 'VVEL_xmean'
    
    SALT_prime = (SALT - SALT_xmean)  # azonal salt component (x,y,z) [g/kg]
    VVEL_prime = (VVEL - VVEL_xmean)  # azonal velocity comp. (x,y,z) [cm/s]
    
    # TRANSPORT terms
    # can integrate now vertically with dz, as transport terms are reduced by PBC_frac
    Fov = ( -1/S0*(VVEL_xint*(SALT_xmean - S0)*ds.dz).sum(dim='z_t'))/1e12  # 1 Sv = 1e12 cm^3/s
    Sov = (VVEL_xint*SALT_xmean*ds.dz).sum(dim='z_t')*rho_sw/1e9  # 1 kg/s = rho_w * 1e-9 g/kg cm^3/s
    
    vSp = (ds.DXU*ds.dz*VVEL_prime*SALT_prime).sum(dim=['nlon_u','z_t'])  # product of primed velocity and salinity [cm^3/s * g/kg]
    Saz = vSp*rho_sw/1e9
    Faz = -vSp/1e12/S0
    
    # total salt transport on (nlon_t, nlat_u)-gridpoints
    # (see /home/ajueling/CESM/doc/how_VNS_is_calculated.txt)
    rn, ac = {'nlat_t':'nlat_u'}, {'nlat_u':ds['nlat_u'].values}
    TVOL_tu = (ds.DZT*ds.TAREA).rename(rn).assign_coords(ac)  # volume T-cell shifted northward to (nlon_u, nlat_t)
    St = rho_sw/1e9*(ds.VNS*TVOL_tu).where(grid.interp(MASK_tt,'Y')).sum(dim=['z_t','nlon_t'])  
    Se = St - rho_sw/1e9*(ds.VVEL*SALT*ds.DXU*ds.DZU).sum(dim=['nlon_u','z_t'])
    
    Fov.name, Faz.name, Sov.name, Saz.name, Se.name, St.name = 'Fov', 'Faz', 'Sov', 'Saz', 'Se', 'St'
    ds = xr.merge([SALT_xmean, VVEL_xmean, Fov, Faz, Sov, Saz, Se, St])
    return ds


def FW_SALT_flux_dataset(run):
    """ combines annual files inot single xr dataset """
    if run=='ctrl':  yy = np.arange(1,301)
    elif run=='lpd':  yy = np.arange(154,601)
    elif run in ['rcp', 'lr1']:  yy = np.arange(2000,2101)

    fn = f'{path_prace}/Mov/FW_SALT_fluxes_{run}.nc'
    if os.path.exists(fn):
        ds = xr.open_dataset(fn, decode_times=False)
    else:
        ds_list = []
        for y in yy:
            ds_ = xr.open_dataset(f'{path_prace}/Mov/FW_SALT_fluxes_{run}_{y:04d}.nc')
            ds_list.append(ds_.copy())
        ds = xr.concat(ds_list, dim='time')
        ds.to_netcdf(fn)
    return ds


if __name__=='__main__':
    """
    input:
    run .. file name

    output:
    ds .. dataset containing FW/SALt transport terms
          saved to f'{path_prace}/Mov/FW_SALT_fluxes_{run}_{y}.nc'
    """
    print(datetime.now())

    run = sys.argv[1]
    ys, ye = int(sys.argv[2]), int(sys.argv[3])

    if run in ['ctrl', 'rcp']:
        fe = file_ex_ocn_ctrl
        domain = 'ocn'
    elif run in ['lpd', 'lr1']:
        fe = file_ex_ocn_lpd
        domain = 'ocn_low'
    geometrics = ['TAREA', 'UAREA', 'dz', 'DXT', 'DXU', 'HTN', 'HUS', 'DYT', 'DYU', 'HTE', 'HUW', 'REGION_MASK']
    ds_geo = xr.open_dataset(fe, decode_times=False)[geometrics].drop(['TLONG','TLAT','ULONG','ULAT']).squeeze()
    DZT = xr_DZ_xgcm(domain=domain, grid='T')
    DZU = xr_DZ_xgcm(domain=domain, grid='U')

    for y in tqdm.tqdm(np.arange(ys,ye)):
        fn = f'{path_prace}/Mov/FW_SALT_fluxes_{run}_{y:04d}.nc'
        if os.path.exists(fn):  print(y, fn, ' exists'); continue
        print(y, fn)
        da_SALT = xr.open_dataset(f'{path_prace}/{run}/ocn_yrly_SALT_{y:04d}.nc', decode_times=False)
        da_VNS  = xr.open_dataset(f'{path_prace}/{run}/ocn_yrly_VNS_{y:04d}.nc' , decode_times=False)
        da_VVEL = xr.open_dataset(f'{path_prace}/{run}/ocn_yrly_UVEL_VVEL_{y:04d}.nc', decode_times=False).VVEL
        
        # some files were created without setting `decode_times` to `False`
        # this creates a `time` variable which messes with xr.merge
        if 'time' in da_SALT:  da_SALT = da_SALT.drop('time') 
        if 'time' in da_VNS:  da_VNS = da_VNS.drop('time') 
        
        ds = xr.merge([da_SALT, da_VNS, da_VVEL, ds_geo])
        ds.VNS .attrs['grid_loc'] = 3121
        ds.VVEL.attrs['grid_loc'] = 3221
        ds.SALT.attrs['grid_loc'] = 3111

        # prepare grid
        (g_, ds_) = pop_tools.to_xgcm_grid_dataset(ds)
        ds_['DZT'] = DZT
        ds_['DZU'] = DZU
        grid = xgcm.Grid(ds_, metrics=metrics, coords=coords)

        MASK_tt = ds_.REGION_MASK
        Atl_MASK_tt = xr.DataArray(np.in1d(MASK_tt, [6,8,9]).reshape(MASK_tt.shape),
                                dims=MASK_tt.dims, coords=MASK_tt.coords)
        rn = {'nlat_t':'nlat_u', 'nlon_t':'nlon_u'}
        ac = {'nlat_u':ds_['nlat_u'].values, 'nlon_u':ds_['nlon_u'].values}
        Atl_MASK_uu = Atl_MASK_tt.rename(rn).assign_coords(ac)

        # calculate transport terms
        ds = calc_transports_pbc(grid=grid, ds=ds_, MASK_tt=Atl_MASK_tt, MASK_uu=Atl_MASK_uu)
        ds.to_netcdf(fn)

    print(datetime.now())
    