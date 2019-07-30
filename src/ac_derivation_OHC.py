import os
import numpy as np
import xarray as xr
import datetime

from paths import path_samoc
from regions import boolean_mask, regions_dict
from constants import cp_sw, rho_sw, km
from timeseries import IterateOutputCESM
from xr_integrate import xr_int_global, xr_int_global_level,\
                         xr_int_zonal, xr_int_zonal_level
from xr_DataArrays import xr_DZ, xr_AREA, xr_HTN, xr_LATS

class DeriveOHC(object):
    """ generate (detrended) OHC files """
    def __init__(self):
        return
    
    
    def generate_OHC_files(self, run, parallel=0):
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

        if run in ['ctrl', 'rcp']:   domain = 'ocn'
        elif run in ['lpd', 'lpi']:  domain = 'ocn_low'

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
        if run in ['ctrl', 'rcp']:
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
            if parallel>0:
                if (y+parallel)%10==0: pass
                else: continue
            
            file_out = f'{path_samoc}/OHC/OHC_integrals_{run}_{y}.nc'

#             if os.path.exists(file_out):
#     #             should check here if all the fields exist
#                 print(f'{datetime.datetime.now()} {y} skipped as files exists already')
#                 continue
            print(f'{datetime.datetime.now()} {y}, {file}')

            t   = y*365  # time in days since year 0, for consistency with CESM date output
            ds  = xr.open_dataset(file, decode_times=False).drop(['ULONG', 'ULAT'])
            round_tlatlon(ds)

            if ds.PD[0,150,200].round(decimals=0)==0:
                ds['PD'] = ds['PD']*1000 + rho_sw
            elif ds.PD[0,150,200].round(decimals=0)==1:
                ds['PD'] = ds['PD']*1000
            else: 
                print('density [g/cm^3] is neither close to 0 or 1')

            OHC = ds.TEMP*ds.PD*cp_sw
#             OHC = ds.TEMP*rho_sw*cp_sw
            ds.close()
            OHC = OHC.where(MASK)

            OHC_DZT = OHC*DZT
            print(f'{datetime.datetime.now()}  {y} calculated OHC & OHC_DZT')
            
            # global, global levels, zonal, zonal levels integrals for different regions
            for mask_nr in [0,1,2,3,6,10]:
                print(datetime.datetime.now(), mask_nr)
                da = OHC.where(boolean_mask(domain, mask_nr=mask_nr))
                print(datetime.datetime.now(), mask_nr)
                
                da_g  = xr_int_global(da=da, AREA=AREA, DZ=DZT)
                print(datetime.datetime.now(), mask_nr)
                da_gl = xr_int_global_level(da=da, AREA=AREA, DZ=DZT)
                print(datetime.datetime.now(), mask_nr)
                da_z  = xr_int_zonal(da=da, HTN=HTN, LATS=LATS, AREA=AREA, DZ=DZT)
                print(datetime.datetime.now(), mask_nr)
                da_zl = xr_int_zonal_level(da=da, HTN=HTN, LATS=LATS, AREA=AREA, DZ=DZT)
                print(datetime.datetime.now(), mask_nr)
                
                name = regions_dict[mask_nr]
                ds_g  = t2ds(da_g , f'OHC_{name}'              , t)
                ds_gl = t2ds(da_gl, f'OHC_levels_{name}'       , t)
                ds_z  = t2ds(da_z , f'OHC_zonal_{name}'        , t)
                ds_zl = t2ds(da_zl, f'OHC_zonal_levels_{name}' , t)
                
                if mask_nr==0:
                    ds_new = xr.merge([ds_g, ds_gl, ds_z, ds_zl])
                else:
                    ds_new = xr.merge([ds_new, ds_g, ds_gl, ds_z, ds_zl])
            print(f'{datetime.datetime.now()}  done with calculations')
            
            # vertical integrals
            da_v  = OHC_DZT.sum(dim='z_t')                         #   0-6000 m
            da_va = OHC_DZT.isel(z_t=slice( 0, 9)).sum(dim='z_t')  #   0- 100 m
            da_vb = OHC_DZT.isel(z_t=slice( 9,20)).sum(dim='z_t')  # 100- 700 m
            da_vc = OHC_DZT.isel(z_t=slice(20,26)).sum(dim='z_t')  # 700-2000 m
            
            ds_v  = t2ds(da_v , 'OHC_vertical_0_6000m'  , t)
            ds_va = t2ds(da_va, 'OHC_vertical_0_100m'   , t)
            ds_vb = t2ds(da_vb, 'OHC_vertical_100_700m' , t)
            ds_vc = t2ds(da_vc, 'OHC_vertical_700_2000m', t)

            ds_new = xr.merge([ds_new, ds_v, ds_va, ds_vb, ds_vc])

            print(f'{datetime.datetime.now()}  done making datasets')
            print(f'output: {file_out}\n')
            ds_new.to_netcdf(path=file_out, mode='w')
            ds_new.close()

            if y in [2002, 102, 156, 1602]:  break  # for testing only

        # combining yearly files
        if parallel==0:
            file_out = f'{path_samoc}/OHC/OHC_integrals_{run}.nc'
            mfname = f'{path_samoc}/OHC/OHC_integrals_{run}_*.nc'
            if os.path.isfile(file_out):  os.remove(file_out)
            combined = xr.open_mfdataset(mfname,
                                         concat_dim='time',
        #                                  autoclose=True,
                                         coords='minimal')
            combined.to_netcdf(file_out)
    #         os.remove(mfname)
        print(f'{datetime.datetime.now()}  done\n')
        
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