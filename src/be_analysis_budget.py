import os
import dask
import numpy as np
import xarray as xr

from tqdm import tqdm
from paths import path_prace
from regions import regions_dict, boolean_mask
from constants import rho_sw, cp_sw
from xr_DataArrays import xr_DZ, xr_AREA, dll_dims_names
from timeseries import IterateOutputCESM


neighbours = [(1,2), (1,3), (1,6),
              (2,3), (2,10),
              (3,5),
              (6,7), (6,8), (6,9),
              (8,10), (8,11),
              (9,10),
             ]


class AnalyzeBudget(object):
    """ calculating heat fluxes """
    
    def __init__(self):
        return
    
    
    def make_adv_cell_dataset(self, run):
        """ dataset that contains advection cells """
#         for pair in neighbours:
            
        return
        
    
    def advection_cells(self, from_basin_mask, to_basin_mask):
        """ arrays with which east-/northward advection need to be multiplied 
        adv_E:   1   if to_basin to the East of from_basin
                -1   if to_basin to the West of from_basin
               nan   elsewhere
        adv_N:   1   if to_basin to the North of from_basin
                -1   if to_basin to the South of from_basin
               nan   elsewhere

        """
        assert np.shape(from_basin_mask)==np.shape(to_basin_mask)
        assert from_basin_mask.dims==to_basin_mask.dims
        (lat, lon) = from_basin_mask.dims
        m0, m1 = from_basin_mask, to_basin_mask
        adv_E = (+m0*m1.roll(shifts={lon:-1}, roll_coords=lon)\
                 -m1*m0.roll(shifts={lon:-1}, roll_coords=lon))\
                .fillna(0)
        adv_N = (+m0*m1.shift(shifts={lat:-1})\
                 -m1*m0.shift(shifts={lat:-1}))\
                .fillna(0)
        if np.all(np.isnan(adv_E)) and np.all(np.isnan(adv_N)):
            print('warning, no neighbouring cells!')
        return adv_E, adv_N
    
    
    def all_advection_cells(self, domain):
        fn = f'{path_prace}/OHC/advection_cells_{domain}.nc'
        try:
            assert os.path.exists(fn), f'{fn} does not exist'
            adv = xr.open_dataset(fn)
        except:
            for i, pair in enumerate(tqdm(neighbours)):
                from_basin_mask = boolean_mask(domain, pair[0])
                to_basin_mask   = boolean_mask(domain, pair[1])
                name = f'{regions_dict[pair[0]]}_to_{regions_dict[pair[1]]}'
                adv_E, adv_N = self.advection_cells(from_basin_mask, to_basin_mask)
                adv_E.name = f'adv_E_{name}'
                adv_N.name = f'adv_N_{name}'
                if i==0: adv = xr.merge([adv_N, adv_E])
                else: adv = xr.merge([adv, adv_N, adv_E])
            adv.to_netcdf(fn)
        return adv
    
    
    def all_transports(self, run, quantity):
        """ computes heat or salt fluxes """
        assert run in ['ctrl', 'lpd']
        

        assert quantity in ['SALT', 'OHC']

        if quantity=='OHC':
            VN, UE = 'VNT', 'UET'
            conversion = rho_sw*cp_sw
            qstr = 'heat'
            unit_out = 'W'
        elif quantity=='SALT':
            VN, UE = 'VNS', 'UES'
            conversion = rho_sw*1e-3
            qstr = 'salt'
            unit_out = 'kg/s'
            
        if run=='ctrl':
            domain='ocn'
            all_transports_list = []
            
        elif run=='lpd':
            domain='ocn_low'
            mf_fn = f'{path_prace}/{run}/ocn_yrly_{VN}_{UE}_*.nc'
            kwargs = {'concat_dim':'time',
                      'decode_times':False,
                      'drop_variables':['TLONG', 'TLAT', 'ULONG', 'ULAT'],
                      'parallel':True}
            ds = xr.open_mfdataset(mf_fn, **kwargs)

        DZ = xr_DZ(domain=domain)
        adv = self.all_advection_cells(domain=domain)
        AREA = xr_AREA(domain=domain)
        dims = [dim for dim in dll_dims_names(domain=domain)]

        for i, pair in enumerate(tqdm(neighbours)):
            name = f'{qstr}_flux_{regions_dict[pair[0]]}_to_{regions_dict[pair[1]]}'
#             if i>2:  continue
            adv_E = adv[f'adv_E_{regions_dict[pair[0]]}_to_{regions_dict[pair[1]]}']
            adv_N = adv[f'adv_N_{regions_dict[pair[0]]}_to_{regions_dict[pair[1]]}']
            MASK  = ((abs(adv_E)+abs(adv_N))/(abs(adv_E)+abs(adv_N))).copy()
            adv_E = adv_E.where(MASK==1, drop=True)
            adv_N = adv_N.where(MASK==1, drop=True)
            DZ_    = DZ.where(MASK==1, drop=True)
            AREA_  = AREA.where(MASK==1, drop=True)
            if run=='ctrl':
                for j, (y,m,f) in tqdm(enumerate(IterateOutputCESM(domain='ocn', run='ctrl',\
                                                                   tavg='yrly', name=f'{VN}_{UE}'))):
                    if j>1: continue
                    ds = xr.open_dataset(f, decode_times=False).where(MASK==1, drop=True)
                    transport = ((adv_E*ds[UE] + adv_N*ds[VN])*AREA_*DZ_).sum(dim=dims)*conversion
                    transport.name = name
                    transport.attrs['units'] = unit_out
                    if j==0:  transport_t = transport
                    else:     transport_t = xr.concat([transport_t, transport], dim='time')
                
                all_transports_list.append(transport_t)
                    
            elif run=='lpd':
                ds_ = ds.where(MASK==1, drop=True)
                transport = ((adv_E*ds_[UE] + adv_N*ds_[VN])*AREA_*DZ_).sum(dim=dims)*conversion
                
                transport.name = name
                transport.attrs['units'] = unit_out
                if i==0: all_transports = transport
                else: all_transports = xr.merge([all_transports, transport])
        
        if run=='ctrl':  all_transports = xr.merge(all_transports_list)
                
        all_transports.to_netcdf(f'{path_prace}/{quantity}/{quantity}_fluxes_{run}.nc')
        return all_transports
    
    
    def surface_heat_flux(self, run):
        """ total surface heat flux into ocean basins """
        # 32:20 min ctrl
        # 1min 4s lpd
        if run=='ctrl':   domain = 'ocn'
        elif run=='lpd':  domain = 'ocn_low'
        da = xr.open_mfdataset(f'{path_prace}/{run}/ocn_yrly_SHF_0*.nc', concat_dim='time').SHF
        AREA = xr_AREA(domain=domain)
        SHF = spy*(da*AREA).sum(dim=['nlat', 'nlon'])
        SHF.name = 'Global_Ocean'
        for nr in tqdm(np.arange(1,12)):
            MASK = boolean_mask(domain=domain, mask_nr=nr)
            temp = spy*(da*AREA).where(MASK).sum(dim=['nlat', 'nlon'])
            temp.name = regions_dict[nr]
            SHF = xr.merge([SHF, temp])
        SHF.attrs['quantity'] = 'yrly averaged total surface heat flux, positive down'
        SHF.attrs['units'] = '[J/yr]'
        SHF.to_netcdf(f'{path_prace}/OHC/SHF_{run}.nc')
        return