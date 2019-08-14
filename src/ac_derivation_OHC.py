import os
import numpy as np
import xarray as xr
import datetime

from tqdm import tqdm
from paths import path_samoc
from regions import boolean_mask, regions_dict
from constants import cp_sw, rho_sw, km
from timeseries import IterateOutputCESM
from xr_integrate import xr_int_global, xr_int_global_level,\
                         xr_int_zonal, xr_int_zonal_level
from xr_DataArrays import xr_DZ, xr_AREA, xr_HTN, xr_LATS, dll_dims_names
from xr_regression import xr_quadtrend

import matplotlib.pyplot as plt

class DeriveOHC(object):
    """ generate (detrended) OHC files """
    def __init__(self):
        return
    
    
    def generate_OHC_files(self, run, year=None):
        """ non-detrended OHC files for full length of simulations
        
        One file contains integrals (all global and by basin):
        
        x,y,z .. scalars        
        x,y   .. vertical profiles 
        x     .. "zonal" integrals 
        
        A separate file each for 4 different depth levels
        z     .. 2D maps, global only, but for different vertical levels
 
        # (ocn:      takes about 45 seconds per year: 70 yrs approx 55 mins)
        (ocn:      takes about 14 min per year)
        (ocn_rect: takes about  3 seconds per year: 70 yrs approx 3 mins)
        """
        
        def t2da(da, t):
            """adds time dimension to xr DataArray, then sets time value to t"""
            da = da.expand_dims('time')
            da = da.assign_coords(time=[t])
            return da

        def t2ds(da, name, t):
            """ 
            adds time dimension to xr DataArray, then sets time value to t,
            and then returns as array in xr dataset
            """
            da = t2da(da, t)
            ds = da.to_dataset(name=name)
            return ds
        
        print(f'\n{datetime.datetime.now()}  start OHC calculation: run={run}')
        assert run in ['ctrl', 'rcp', 'lpd', 'lpi']

        if run=='rcp':
            domain = 'ocn'
        elif run=='ctrl':
            domain = 'ocn_rect'
        elif run in ['lpd', 'lpi']:
            domain = 'ocn_low'
            
        (z, lat, lon) = dll_dims_names(domain)

        # geometry
        DZT  = xr_DZ(domain)
        AREA = xr_AREA(domain)
        HTN  = xr_HTN(domain)
        LATS = xr_LATS(domain)
        
        def round_tlatlon(das):
            """ rounds TLAT and TLONG to 2 decimals
            some files' coordinates differ in their last digit
            rounding them avoids problems in concatonating
            """
            das['TLAT']   = das['TLAT'].round(decimals=2)
            das['TLONG']  = das['TLONG'].round(decimals=2)
            return das
        if domain=='ocn':
            round_tlatlon(HTN)
            round_tlatlon(LATS)

        MASK = boolean_mask(domain, mask_nr=0)
        for k in range(42):
            DZT[k,:,:]  = DZT[k,:,:].where(MASK)
        AREA = AREA.where(MASK)
        HTN  = HTN.where(MASK)
        LATS = LATS.where(MASK)
        print(f'{datetime.datetime.now()}  done with geometry')

        for y,m,file in IterateOutputCESM(domain, run, 'yrly', name='TEMP_PD'):
            
            if year!=None:  # select specific year
                if year==y:
                    pass
                else:
                    continue
            
            file_out = f'{path_samoc}/OHC/OHC_integrals_{run}_{y:04d}.nc'

#             if os.path.exists(file_out):
#     #             should check here if all the fields exist
#                 print(f'{datetime.datetime.now()} {y} skipped as files exists already')
#                 continue
            print(f'{datetime.datetime.now()} {y}, {file}')

            t   = y*365  # time in days since year 0, for consistency with CESM date output
            ds  = xr.open_dataset(file, decode_times=False).TEMP
            if domain=='ocn':
                ds = ds.drop(['ULONG', 'ULAT'])
                ds = round_tlatlon(ds)

#             if ds.PD[0,150,200].round(decimals=0)==0:
#                 ds['PD'] = ds['PD']*1000 + rho_sw
#             elif ds.PD[0,150,200].round(decimals=0)==1:
#                 ds['PD'] = ds['PD']*1000
#             else: 
#                 print('density [g/cm^3] is neither close to 0 or 1')

#             OHC = ds.TEMP*ds.PD*cp_sw
            OHC = ds*rho_sw*cp_sw
            ds.close()
            OHC = OHC.where(MASK)

            OHC_DZT = OHC*DZT
            print(f'{datetime.datetime.now()}  {y} calculated OHC & OHC_DZT')
            
            # global, global levels, zonal, zonal levels integrals for different regions
            for mask_nr in tqdm([0,1,2,3,6,7,8,9,10]):
                name = regions_dict[mask_nr]
                da = OHC.where(boolean_mask(domain, mask_nr=mask_nr))
                
                da_g = (da*AREA*DZT).sum(dim=[z, lat, lon])
                da_g.attrs['units'] = '[J]'
                ds_g  = t2ds(da_g , f'OHC_{name}', t)

                da_gl = (da*AREA).sum(dim=[lat, lon])
                da_gl.attrs['units'] = '[J m^-1]'
                ds_gl = t2ds(da_gl, f'OHC_levels_{name}', t)

                if domain=='ocn':  da_z  = xr_int_zonal(da=da, HTN=HTN, LATS=LATS, AREA=AREA, DZ=DZT)
                else:  da_z = (da*HTN*DZT).sum(dim=[z, lon])
                da_z.attrs['units'] = '[J m^-1]'
                ds_z = t2ds(da_z , f'OHC_zonal_{name}', t)
                
                if domain=='ocn':  da_zl = xr_int_zonal_level(da=da, HTN=HTN, LATS=LATS, AREA=AREA, DZ=DZT)
                else:  da_zl = (da*HTN).sum(dim=[lon])
                da_zl.attrs['units'] = '[J m^-2]'
                ds_zl = t2ds(da_zl, f'OHC_zonal_levels_{name}', t)
                if mask_nr==0:   ds_new = xr.merge([ds_g, ds_gl, ds_z, ds_zl])
                else:            ds_new = xr.merge([ds_new, ds_g, ds_gl, ds_z, ds_zl])
                    
            print(f'{datetime.datetime.now()}  done with horizontal calculations')
            
            # vertical integrals
            # full depth
            da_v  = OHC_DZT.sum(dim=z)                         #   0-6000 m
            da_v.attrs = {'depths':f'{OHC_DZT[z][0]-OHC_DZT[z][-1]}',
                          'units':'[J m^-2]'}
            
            if domain in ['ocn', 'ocn_rect']:  zsel = [[0,9], [0,20], [20,26]]
            elif domain=='ocn_low':            zsel = [[0,9], [0,36], [36,45]]
            
            #   0- 100 m
            da_va = OHC_DZT.isel({z:slice(zsel[0][0], zsel[0][1])}).sum(dim=z)  
            da_va.attrs = {'depths':f'{OHC_DZT[z][zsel[0][0]].values:.0f}-{OHC_DZT[z][zsel[0][1]].values:.0f}',
                           'units':'[J m^-2]'}
            
            #   0- 700 m
            da_vb = OHC_DZT.isel({z:slice(zsel[1][0],zsel[1][1])}).sum(dim=z)  
            da_vb.attrs = {'depths':f'{OHC_DZT[z][zsel[1][0]].values:.0f}-{OHC_DZT[z][zsel[1][1]].values:.0f}',
                           'units':'[J m^-2]'}
            
            # 700-2000 m
            da_vc = OHC_DZT.isel({z:slice(zsel[2][0],zsel[2][1])}).sum(dim=z)  
            da_vc.attrs = {'depths':f'{OHC_DZT[z][zsel[2][0]].values:.0f}-{OHC_DZT[z][zsel[2][1]].values:.0f}',
                           'units':'[J m^-2]'}
            
            ds_v  = t2ds(da_v , 'OHC_vertical_0_6000m'  , t)
            ds_va = t2ds(da_va, 'OHC_vertical_0_100m'   , t)
            ds_vb = t2ds(da_vb, 'OHC_vertical_0_700m'   , t)
            ds_vc = t2ds(da_vc, 'OHC_vertical_700_2000m', t)

            ds_new = xr.merge([ds_new, ds_v, ds_va, ds_vb, ds_vc])

            print(f'{datetime.datetime.now()}  done making datasets')
            print(f'output: {file_out}\n')
            ds_new.to_netcdf(path=file_out, mode='w')
            ds_new.close()

#             if y in [2002, 102, 156, 1602]:  break  # for testing only

        # combining yearly files
        
        print(f'{datetime.datetime.now()}  done\n')
        
        if run=='ctrl':  print('year 205 is wrong and should be averaged by executing `fix_ctrl_year_205()`')
        return
    
    
    def combine_yrly_OHC_integral_files(self, run):
        file_out = f'{path_samoc}/OHC/OHC_integrals_{run}.nc'
        mfname = f'{path_samoc}/OHC/OHC_integrals_{run}_*.nc'
#         if os.path.isfile(file_out):  os.remove(file_out)
        combined = xr.open_mfdataset(mfname,
                                     concat_dim='time',
        #                                  autoclose=True,
                                     coords='minimal')
        combined.to_netcdf(file_out)
#         if os.path.isfile(file_out): os.remove(mfname)
        return
    
    
    def fix_ctrl_year_205():
        """ averages year 204 and 206 to replace year 205 in ctrl run """
        ctrl = xr.open_dataset(f'{path_samoc}/OHC/OHC_integrals_ctrl.nc')
        for field in ['', '_levels', '_zonal', '_zonal_levels']:
            for ocean in ['Global', 'Atlantic', 'Southern', 'Indian', 'Pacific', 'Arctic']:
                ctrl[f'OHC{field}_{ocean}_Ocean'][204] = (ctrl[f'OHC{field}_{ocean}_Ocean'][203] + ctrl[f'OHC{field}_{ocean}_Ocean'][205])/2
        for depths in ['0_6000m', '0_100m', '0_700m', '700_2000m']:
            ctrl[f'OHC_vertical_{depths}'][204] = (ctrl[f'OHC_vertical_{depths}'][204] + ctrl[f'OHC_vertical_{depths}'][204])/2
        ctrl.to_netcdf(f'{path_samoc}/OHC/OHC_integrals_ctrl_.nc')
        ctrl.close()
        cmd = f'mv {path_samoc}/OHC/OHC_integrals_ctrl_.nc {path_samoc}/OHC/OHC_integrals_ctrl.nc '
        os.system(cmd)
        return
    
    def quadratically_detrend_OHC_integral_fields(ds, name):
        """ iterates through all fields in OHC_integral file
        and quadratically detrends them in time
        """
        # < 1 min for both ctrl_rect and lpd
        ds_ = ds.copy()
        key_list = [k for k in ds.data_vars.keys()]
        for k in tqdm(key_list):
            ds_[k].values = (ds[k] - xr_quadtrend(ds[k])).values
        ds_.to_netcdf(f'{path_samoc}/OHC/OHC_integrals_{name}_qd.nc')
        return
        
        
    def OHC_pointwise_detrending(self, run):
        """ removes quadratic trend at each point """
        return
        
        
if __name__=="__main__":
    run      = sys.argv[1]
    parallel = int(sys.argv[2])

    assert run in ['ctrl', 'rcp', 'lpd', 'lpi']
    
    # 15 min for one lpi run
  
#     if run in ['lpd', 'lpi']:     # low res
    
#         cluster = LocalCluster(n_workers=4)
#         client = Client(cluster)
#         OHC_parallel(run=run, mask_nr=mask_nr)
        
#     elif run in ['ctrl', 'rcp']:  # high res
#         OHC_integrals(run=run, mask_nr=mask_nr)
    DeriveOHC().generate_OHC_files(run=run, parallel=parallel)
    
    print(f'\n\nfinished at\n{datetime.datetime.now()}')