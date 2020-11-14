""" Freshwater/Salt Transport
- freshwater relative to 35 g/kg
- taking partial bottom cells of high resolution POP grid into account

if called directly, this scripts calculates:

(1) the total salt/volume/freshwater flux through Bering Strait and Strait of Gibraltar
->  f'{path_prace}/Mov/FW_Med_Bering_{run}_{y:04d}.nc'

(2) the Atlantic freshwater and salt meridional transport terms
    (overturning, azonal, total, and eddy) as a function of latitude
->  f'{path_prace}/Mov/FW_SALT_fluxes_{run}_{y:04d}.nc'
"""
import os
import sys
import tqdm
import xgcm
import dask
import numpy as np
import xarray as xr
import warnings
import pop_tools

warnings.filterwarnings('ignore')

from paths import file_ex_ocn_ctrl, file_ex_ocn_lpd, path_prace, path_results
from datetime import datetime
from constants import rho_sw
from xr_DataArrays import xr_DZ, xr_DZ_xgcm

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
sections_high = {
    'Bering': dict(nlon=slice(2929,2947), nlat=1979),
    'Med'   : dict(nlon=1042, nlat=slice(1555,1559))
}
sections_low = {
    'Bering': dict(nlon=slice(198 ,202 ), nlat=333),
    'Med':    dict(nlon=32,nlat=slice(290,295))
}


def calc_section_transport(ds):
    """ calculate the transport across Bering Strait and Gibraltar sections
    section  predermined (nlon,nlat) dictionary
    """
    for q in ['VVEL', 'UVEL', 'VNS', 'UES', 'DXU', 'DYU', 'TAREA', 'DZT', 'DZU']:
        assert q in ds
        
    if ds.TAREA.shape==(2400, 3600):  section = sections_high
    elif ds.TAREA.shape==(384, 320):  section = sections_low
    
    S0 = 35.  # [g/kg]
    
    S_BS = (ds.VNS*ds.TAREA*ds.DZT).isel(section['Bering']).sum()/1e4  # [cm^2*m/s] -> [kg_{Salt}/s]
    V_BS = (ds.VVEL*ds.DXU*ds.DZU).isel(section['Bering']).sum()/1e4   # [cm^2*m/s] -> [m^3/s]
    F_BS = -1/S0*(S_BS - S0*V_BS)                                   # [m^3/s]
    
    S_Med = (ds.UES*ds.TAREA*ds.DZT).isel(section['Med']).sum()/1e4    # [cm^2*m/s] -> [kg_{Salt}/s]
    V_Med = (ds.UVEL*ds.DYU*ds.DZU).isel(section['Med']).sum()/1e4     # [cm^2*m/s] -> [m^3/s]
    F_Med = -1/S0*(S_Med - S0*V_Med)

    S_BS.name, V_BS.name, F_BS.name, = 'S_BS', 'V_BS', 'F_BS'
    S_Med.name, V_Med.name, F_Med.name, = 'S_Med', 'V_Med', 'F_Med'
    
    S_BS.attrs['units'], S_Med.attrs['units'] = 'kg/s', 'kg/s'
    V_BS.attrs['units'], V_Med.attrs['units'] = 'm^3/s', 'm^3/s'
    F_BS.attrs['units'], F_Med.attrs['units'] = 'm^3/s', 'm^3/s'
    
    return xr.merge([S_BS, V_BS, F_BS, S_Med, V_Med, F_Med]).expand_dims('time')


def calc_transports_pbc(grid, ds, MASK_tt, MASK_uu, S0=35.0):
    """ calculate the merid. FW and SALT transport terms with xgcm grid including PBCs """
    # SALT
    SALT = grid.interp(grid.interp(ds.SALT.where(MASK_tt).where(ds.SALT>0), 'Y'), 'X')  # SALT interpolated to UU grid
    VVEL = (ds.VVEL*ds.DZU/ds.dz).where(MASK_uu).where(ds.VVEL<1000)  # reduced velocity due to PBCs

    # non-PBC-weighted mean here, as salinity differences between neighbouring cells are assumed very small
    SALT_xmean = SALT.mean('nlon_u')
    SALT_xmean.name = 'SALT_xmean'
    
    # VVEL terms:  weighted by cell thickness
    VVEL_baro  = (VVEL*ds.DXU*ds.DZU).sum(['nlon_u','z_t'])/(ds.DXU*ds.DZU).where(VVEL<1000).sum(['nlon_u','z_t'])  # barotropic velocity
    VVEL_xint  = (VVEL*ds.DXU).sum('nlon_u')  # zonal integral velocity (y,z) [cm^2/s]
    VVEL_xint_s  = ((VVEL-VVEL_baro)*ds.DXU).sum('nlon_u')  # zonal integral of pure overturning velocity (y,z) [cm^2/s]
    VVEL_xmean = (VVEL*ds.DXU).sum('nlon_u')/ds.DXU.where(VVEL<1000).sum('nlon_u')  # zonal mean velocity     (y,z) [cm/s]
    VVEL_baro.name = 'VVEL_barotropic'
    VVEL_xmean.name = 'VVEL_xmean'
    
    SALT_prime = (SALT - SALT_xmean)  # azonal salt component (x,y,z) [g/kg]
    VVEL_prime = (VVEL - VVEL_xmean)  # azonal velocity comp. (x,y,z) [cm/s]
    
    # TRANSPORT terms
    # can integrate now vertically with dz, as transport terms are reduced by PBC_frac
    Fov = ( -1/S0*(VVEL_xint*(SALT_xmean - S0)*ds.dz).sum(dim='z_t'))/1e12  # 1 Sv = 1e12 cm^3/s
    Sov = (VVEL_xint*SALT_xmean*ds.dz).sum(dim='z_t')*rho_sw/1e9  # 1 kg/s = rho_w * 1e-9 g/kg cm^3/s
    
    Fovs = ( -1/S0*(VVEL_xint_s*(SALT_xmean - S0)*ds.dz).sum(dim='z_t'))/1e12  # 1 Sv = 1e12 cm^3/s
    Sovs = (VVEL_xint_s*SALT_xmean*ds.dz).sum(dim='z_t')*rho_sw/1e9  # 1 kg/s = rho_w * 1e-9 g/kg cm^3/s

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
    Fovs.name, Sovs.name = 'Fovs', 'Sovs'
    ds = xr.merge([SALT_xmean, VVEL_baro, VVEL_xmean, Fov, Fovs, Faz, Sov, Sovs, Saz, Se, St])
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
    input: {run} {ys} {ye}
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
    else:
        raise ValueError(f'`run`={run} not implemented')
    geometrics = ['TAREA', 'UAREA', 'dz', 'DXT', 'DXU', 'HTN', 'HUS', 'DYT', 'DYU', 'HTE', 'HUW', 'REGION_MASK']
    ds_geo = xr.open_dataset(fe, decode_times=False)[geometrics].drop(['TLONG','TLAT','ULONG','ULAT']).squeeze()
    DZT = xr_DZ(domain)
    DZU = xr_DZ(domain, grid='U')
    DZT.name, DZU.name = 'DZT', 'DZU'
    DZTx = xr_DZ_xgcm(domain=domain, grid='T')
    DZUx = xr_DZ_xgcm(domain=domain, grid='U')

    for y in tqdm.tqdm(np.arange(ys,ye)):
        #region: section transport
        fn = f'{path_prace}/Mov/FW_Med_Bering_{run}_{y:04d}.nc'
        if  os.path.exists(fn): #
            # continue 
            pass
            # print(y, fn, ' exists')
        else:
            da_VNS  = xr.open_dataset(f'{path_prace}/{run}/ocn_yrly_VNS_{y:04d}.nc' , decode_times=False)
            da_UES  = xr.open_dataset(f'{path_prace}/{run}/ocn_yrly_UES_{y:04d}.nc' , decode_times=False)
            ds_VEL  = xr.open_dataset(f'{path_prace}/{run}/ocn_yrly_UVEL_VVEL_{y:04d}.nc', decode_times=False)
            
            # some files were created without setting `decode_times` to `False`
            # this creates a `time` variable which messes with xr.merge
            if 'time' in da_VNS:   da_VNS  = da_VNS .drop('time')
            if 'time' in da_UES:   da_UES  = da_UES .drop('time')

            # some UES files have no name
            # if run=='ctrl':  
            #     da_UES  = xr.open_dataarray(f'{path_prace}/{run}/ocn_yrly_UES_{y:04d}.nc' , decode_times=False)
            # if e:  
            #     da_UES  = xr.open_dataarray(f'{path_prace}/{run}/ocn_yrly_UES_{y:04d}.nc' , decode_times=False)
            
            # print("name:   ", da_UES.name)
            ds = xr.merge([da_VNS, da_UES, ds_VEL, ds_geo, DZT, DZU])
            ds = calc_section_transport(ds)
            ds.to_netcdf(fn)
        #endregion

        #region: meridional transport terms
        fn = f'{path_prace}/Mov/FW_SALT_fluxes_{run}_{y:04d}.nc'
        if 1==0:#os.path.exists(fn):
            pass
            # print(y, fn, ' exists')
            
        else:
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
            ds_['DZT'] = DZTx
            ds_['DZU'] = DZUx
            grid = xgcm.Grid(ds_, metrics=metrics, coords=coords)

            MASK_tt = ds_.REGION_MASK
            Atl_MASK_tt = xr.DataArray(np.in1d(MASK_tt, [6,8,9,12]).reshape(MASK_tt.shape),
                                    dims=MASK_tt.dims, coords=MASK_tt.coords)
            rn = {'nlat_t':'nlat_u', 'nlon_t':'nlon_u'}
            ac = {'nlat_u':ds_['nlat_u'].values, 'nlon_u':ds_['nlon_u'].values}
            Atl_MASK_uu = Atl_MASK_tt.rename(rn).assign_coords(ac)

            # calculate transport terms
            ds = calc_transports_pbc(grid=grid, ds=ds_, MASK_tt=Atl_MASK_tt, MASK_uu=Atl_MASK_uu)
            ds.to_netcdf(fn)
        #endregion

    print(datetime.now())
    