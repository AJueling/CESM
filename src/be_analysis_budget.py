import numpy as np
import xarray as xr

from tqdm import tqdm
from paths import path_prace
from regions import regions_dict, boolean_mask


neighbours = [(1,2), (1,3), (1,6),
              (2,3), (2,10),
              (3,5),
              (6,7), (6,8), (6,9),
              (8,10), (8,11)
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
            assert os.path.exists(fn)
            adv = xr.open_dataset(fn)
        except:
            for i, pair in tqdm(enumerate(neighbours)):
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
    
    
    def all_transports(self, domain, basin, VN_adv, UE_adv):
        """ computes heat or salt fluxes """
        print(f'fluxes into {basin}')
        assert domain in ['ocn', 'ocn_low']
        assert VN_adv.units==UE_adv.units
        if VN_adv.units=='degC/s':
            conversion = rho_sw*cp_sw
            unit = 'W'
        elif VN_adv.units=='gram/kilogram/s':
            conversion = rho_sw*1e-3
            unit = 'kg/s'
        else:
            raise ValueError('units need to be in "degC/s" or "gram/kilogram/s"')

        dims = [dim for dim in dll_dims_names(domain=domain)]
        DZ = xr_DZ(domain=domain)
        AREA = xr_AREA(domain=domain)

        adv_E, adv_N = advection_cells(from_basin_mask=neighbour_mask, to_basin_mask=basin_mask)
        transport = ((adv_E*UE_adv + adv_N*VN_adv)*AREA*DZ).sum(dim=dims)*conversion
        transport.name = f'{regions_dict[n]}'
        transport.attrs['units'] = unit
        if i==0: temp=transport
        else: temp = xr.merge([temp, transport])
    
        return temp
    
    
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