import os
import numpy as np
import xarray as xr
import pandas as pd
import statsmodels.api as sm

from tqdm import tqdm

from OHC import t2ds
from paths import CESM_filename, path_samoc, path_prace, path_data, file_ex_ocn_ctrl, file_ex_ocn_lpd
from regions import boolean_mask
from timeseries import IterateOutputCESM
from xr_DataArrays import depth_lat_lon_names, xr_DZ, xr_DXU, xr_AREA
from xr_regression import xr_lintrend, xr_quadtrend
from ba_analysis_dataarrays import AnalyzeDataArray as ADA


# 10x 149 year time segments and one 250 year segment
times_ctrl = (list(zip(np.arange(51, 150, 10), np.arange(200, 301, 10))))
times_lpd  = (list(zip(np.arange(154, 253, 10), np.arange(303, 452, 10))))
times_ctrl.append((51,301))
times_lpd .append((154,404))


class DeriveSST(object):
    """ generate SST fields """
    def __init__(self):
        return

    
    def generate_yrly_SST_files(self, run):
        """ generate the SST data files from TEMP_PD yearly averaged files """
        # ca. 4:30 min for ctrl/rcp, 1:25 for lpi
        # stacking files into one xr DataArray object
        for i, (y,m,s) in enumerate(IterateOutputCESM('ocn', run, 'yrly', name='TEMP_PD')):
            print(y)
            da = xr.open_dataset(s, decode_times=False).TEMP[0,:,:]
            da = da.drop(['z_t', 'ULONG', 'ULAT'])
            da_time = int(da.time.item())
            if run=='ctrl':
                # years 5-50 have different TLAT/TLON grids
                # somehow the non-computed boxes changed (in the continents)
                if i==0:
                    TLAT = da['TLAT'].round(decimals=2)
                    TLONG = da['TLONG'].round(decimals=2)
                da['TLAT'] = TLAT
                da['TLONG'] = TLONG
            else:
                da['TLAT' ] = da['TLAT' ].round(decimals=2)
                da['TLONG'] = da['TLONG'].round(decimals=2)
            del da.encoding["contiguous"]
            ds = t2ds(da=da, name='SST', t=da_time)
            ds.to_netcdf(path=f'{path_samoc}/SST/SST_yrly_{run}_{y:04}.nc', mode='w')
        combined = xr.open_mfdataset(f'{path_samoc}/SST/SST_yrly_{run}_*.nc',
                                     concat_dim='time', autoclose=True, coords='minimal')
        combined.to_netcdf(f'{path_samoc}/SST/SST_yrly_{run}.nc')
        self.remove_superfluous_files(f'{path_samoc}/SST/SST_yrly_{run}_*.nc')
        return
    

    def generate_yrly_SST_ctrl_rect(self):
        """ monthly data of ctrl run """
        monthly_ctrl = xr.open_dataarray(f'{path_prace}/SST/SST_monthly_ctrl.nc')
        t_bins = np.arange(0,len(monthly_ctrl)+1,12)
        yrly_ctrl = monthly_ctrl.groupby_bins('time', t_bins, right=False).mean(dim='time')
        yrly_ctrl = yrly_ctrl.assign_coords(time_bins=np.arange(1, 301)).rename({'time_bins':'time'})
        yrly_ctrl.to_netcdf(f'{path_prace}/SST/SST_yrly_rect_ctrl.nc')
        return
    
    
    def generate_yrly_global_mean_SST(self, run):
        """ calcaultes the global mean sea surface temperature
        ca. 37 sec for ctrl """
        assert run in ['ctrl', 'lpd']
        
        da = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_{run}.nc', decode_times=False)
        if run=='ctrl':
            AREA = xr_AREA(domain='ocn')
            REGION_MASK = xr.open_dataset(file_ex_ocn_ctrl, decode_times=False).REGION_MASK
        elif run=='lpd':
            AREA = xr_AREA(domain='ocn_low')
            REGION_MASK = xr.open_dataset(file_ex_ocn_lpd, decode_times=False).REGION_MASK
        AREA_total = AREA.where(REGION_MASK>0).sum(dim=['nlat', 'nlon'], skipna=True)
        print(AREA_total)
        da_new = (da*AREA).where(REGION_MASK>0).sum(dim=['nlat', 'nlon'], skipna=True)/AREA_total
        da_new.to_netcdf(f'{path_samoc}/SST/GMSST_yrly_{run}.nc')
        return
            
            
    def generate_monthly_SST_files(self, run):
        """ concatonate monthly files, ocn_rect for high res runs"""
        # 8 mins for 200 years of ctrl
        if run in ['ctrl', 'rcp']:  domain = 'ocn_rect'
        elif run in ['lpd', 'lpi']:  domain = 'ocn_low'
            
        for y,m,s in IterateOutputCESM(domain=domain, tavg='monthly', run=run):
            if m==1: print(y)
            if run in ['ctrl', 'rcp']:   xa = xr.open_dataset(s, decode_times=False).TEMP[0,:,:]
            if run in ['lpd', 'lpi']:    xa = xr.open_dataset(s, decode_times=False).TEMP[0,0,:,:]
            if m==1:   xa_out = xa.copy()    
            else:      xa_out = xr.concat([xa_out, xa], dim='time')
            if m==12:  xa_out.to_netcdf(f'{path_prace}/SST/SST_monthly_{run}_{y:04}.nc')
            # this also means only full years are written out
                        
        combined = xr.open_mfdataset(f'{path_prace}/SST/SST_monthly_{run}_*.nc',
                                     concat_dim='time', decode_times=False)
        if run=='ctrl':
            time_ctrl = np.arange(1+1/24, 301, 1/12)
            combined = combined.assign_coords(time=time_ctrl)
            # something is wrong in Febraury 99, so I average 
            n = np.where(np.isclose(time_ctrl, 99.125))[0][0]
            print(f'month 99/2 is averaged; n={n}')
            combined[n] = (combined[n-1]+combined[n+1])/2
        
        combined.to_netcdf(f'{path_prace}/SST/SST_monthly_{run}.nc')
        combined.close()
        self.remove_superfluous_files(f'{path_samoc}/SST/SST_monthly_{run}_*.nc')
        return
    
    
    def deseasonalize_monthly_data(self, run):
        """ """
        assert run in ['ctrl', 'lpd', 'had']
        monthly = xr.open_dataarray(f'{path_prace}/SST/SST_monthly_{run}.nc')
        if run=='ctrl:':   yrly = xr.open_dataarray(f'{path_prace}/SST/SST_yrly_rect_ctrl.nc')
        else:              yrly = xr.open_dataarray(f'{path_prace}/SST/SST_yrly_{run}.nc')
        assert len(monthly)/len(yrly) == 12.0
        
        temp = monthly.copy()
        for j in tqdm(range(12)):
            m = monthly.isel(time=slice(j,len(monthly)+1,12))
            temp[j::12] -= (m-yrly.assign_coords(time=m.time)).mean(dim='time')
        temp.to_netcdf(f'{path_prace}/SST/SST_monthly_deseasonalized_{runs}.nc')
        return
    
    
    def detrend_monthly_data_pointwise(self, run):
        """ quadratically detrend monthly fields """
        assert run in ['ctrl', 'lpd']
        
        da = xr.open_dataarray(f'{path_prace}/SST/SST_monthly_{run}.nc')
        if run=='ctrl':   times = times_ctrl
        elif run=='lpd':  times = times_lpd
            
        for j in tqdm(range(len(times))):
            time = times[j]
        
            if run=='ctrl':     # ctrl
                print(i, j)
                da_sel = da.sel(time=slice(time[0], time[1]))

            elif run=='lpd':   # lpd
                da_sel = da.sel(time=slice(time[0]*365, time[1]*365))#+1e-4))
                if len(da_sel.time) not in [1788, 3000]:
                    da_sel = da.sel(time=slice(time[0]*365, time[1]*365+1e-4))
                y1, y2 = int(time[0]/365), int(time[0]/365)

            print(time, time[0], time[1], len(da_sel.time), da_sel.time[0].values, da_sel.time[-1].values)
            (da_sel-xr_quadtrend(da_sel)).to_netcdf(f'{path_prace}/SST/SST_monthly_ds_dt_{run}_{time[0]}_{time[1]}.nc')
        return
    
        
    def detrend_monthly_obs_two_factor(self):
        """ remove linear combination of anthropogenic and natural forcing signal from CMIP5 MMM """
        MMM_natural = xr.open_dataarray(f'{path_samoc}/GMST/CMIP5_natural.nc', decode_times=False)
        MMM_anthro  = xr.open_dataarray(f'{path_samoc}/GMST/CMIP5_anthro.nc' , decode_times=False)
        monthly_MMM_natural = np.repeat(MMM_natural, 12)
        monthly_MMM_anthro  = np.repeat(MMM_anthro , 12)
        monthly_MMM_natural = monthly_MMM_natural.assign_coords(time=monthly_had.time)
        monthly_MMM_anthro  = monthly_MMM_anthro .assign_coords(time=monthly_had.time)
        
        forcings = monthly_MMM_natural.to_dataframe(name='natural').join(
                    monthly_MMM_anthro.to_dataframe(name='anthro'))

        SST_stacked = monthly_ds_had.stack(z=('latitude', 'longitude'))
        ds_anthro   = SST_stacked[0,:].squeeze().copy()
        ds_natural  = SST_stacked[0,:].squeeze().copy()

        # multiple linear regression
        X = sm.add_constant(forcings[['anthro', 'natural']])
        for i, coordinate in tqdm(enumerate(SST_stacked.z)):
            y = SST_stacked[:, i].values
            model = sm.OLS(y, X).fit()
            ds_anthro[i] = model.params['anthro']
            ds_natural[i] = model.params['natural']

        beta_anthro  = ds_anthro .unstack('z')
        beta_natural = ds_natural.unstack('z')

        ds = xr.merge([{'forcing_anthro': monthly_MMM_anthro}, {'beta_anthro': beta_anthro}])
        ds.to_netcdf(f'{path_prace}/SST/SST_beta_anthro_MMM_monthly_had.nc')

        ds = xr.merge([{'forcing_natural': monthly_MMM_natural}, {'beta_natural':beta_natural}])
        ds.to_netcdf(f'{path_prace}/SST/SST_beta_natural_MMM_monthly_had.nc')
        
        monthly_ds_dt_had = monthly_ds_had.assign_coords(time=monthly_MMM_anthro.time) \
                            - beta_anthro*monthly_MMM_anthro \
                            - beta_natural*monthly_MMM_natural
        monthly_ds_dt_had.to_netcdf(f'{path_prace}/SST/SST_monthly_ds_tfdt_had.nc')
        
        return
        
    
    def isolate_Pacific_SSTs(self, run, extent):
        """"""
        assert run in ['ctrl', 'lpd', 'had']
        assert extent in ['38S', 'Eq', '20N']
    
        if run=='ctrl':
            monthly_fns = [f'{path_prace}/SST/SST_monthly_ds_dt_ctrl_{time[0]}_{time[1]}.nc' for time in times_ctrl]
        elif run=='lpd':
            monthly_fns = [f'{path_prace}/SST/SST_monthly_ds_dt_lpd_{time[0]}_{time[1]}.nc' for time in times_lpd]
        elif run=='had':
            monthly_fns = [f'{path_prace}/SST/SST_monthly_ds_tfdt_had.nc']
            
        area = xr.open_dataarray(f'{path_prace}/geometry/AREA_{extent}_{domain}.nc')  # created in SST_PDO.ipynb
        for k, fn in tqdm(enumerate(monthly_fns)):
            da = xr.open_dataarray(fn)
            if run=='had':  da = self.shift_had(da)
            da = self.focus_data(da)
            da = da.where(area)
            if run in ['ctrl', 'lpd']:
                time = [times_ctrl, times_lpd][j][k]
                fn = f'{path_prace}/SST/SST_monthly_ds_dt_{extent}_{run}_{time[0]}_{time[1]}.nc'
            else:
                fn = f'{path_prace}/SST/SST_monthly_ds_dt_{extent}_{run}.nc'
            da.to_netcdf(fn)
            
        return
    
    
    @staticmethod
    def generate_monthly_regional_SST_files(run):
        """"""
            # 6 mins for 200 years of ctrl
        for i, r in enumerate(['Pac_38S', 'Pac_Eq', 'Pac_20N']):
            latS = [-38, 0, 20][i]
            lonE = [300, 285, 255][i]
            
            if run  in ['ctrl', 'rcp']:  domain= 'ocn_rect'
            elif run in ['lpd', 'lpi']:  domain = 'ocn_low'
                
            Pac_MASK = mask_box_in_region(domain=domain, mask_nr=2,
                                          bounding_lats=(latS,68),
                                          bounding_lons=(110,lonE))
            if run in ['ctrl', 'rcp']:
                Pac_MASK = Pac_MASK.where(Pac_MASK.t_lon+1/.6*Pac_MASK.t_lat<333,0)
            NPac_area = xr_AREA(domain).where(Pac_MASK, drop=True)
            
            for y,m,s in IterateOutputCESM(domain=domain, tavg='monthly', run=run):
                xa = xr.open_dataset(s, decode_times=False).TEMP[0,:,:].where(Pac_MASK, drop=True)
                if m==1:
                    print(y)
                    xa_out = xa.copy()   
                else:
                    xa_out = xr.concat([xa_out, xa], dim='time')
                if m==12:
                    xa_out.to_netcdf(f'{path_samoc}/SST/SST_monthly_{r}_{run}_{y:04}.nc')
                            
            combined = xr.open_mfdataset(f'{path_samoc}/SST/SST_monthly_{r}_{run}_*.nc',
                                         concat_dim='time', decode_times=False)
            combined.to_netcdf(f'{path_samoc}/SST/SST_monthly_{r}_{run}.nc')
            combined.close()
            
#             self.remove_superfluous_files(f'{path_samoc}/SST/SST_monthly_{r}_{run}_*.nc')
            # # remove yearly files
        return
            
        
    @staticmethod
    def generate_monthly_mock_linear_GMST_files(run):
        """ generates xr DataArray with same time coordinates as monthly SST fields
        which contains the linear fit to the mean SST as a stand in for the missing
        monthly GMST
        """
        assert run in ['ctrl', 'lpd']
        if run=='ctrl':   dims = ['t_lon', 't_lat']
        elif run=='lpd':  dims = ['nlon', 'nlat']
        da = xr.open_dataarray(f'{path_samoc}/SST/SST_monthly_{run}.nc', decode_times=False)
        da_new = xr_lintrend(da.mean(dim=dims, skipna=True, keep_attrs=True))
        da_new.name = 'GMST'
        da_new.attrs = {'Note':'This is the linear trend of the SST evolution, not GMST'}
        da_new.to_dataset().to_netcdf(f'{path_samoc}/GMST/GMST_monthly_{run}.nc')
        return
        
        
    def SST_remove_forced_signal(self, run, tres='yrly', detrend_signal='GMST', time_slice='full'):
        """ detrending the SST field
        a) remove the scaled, forced MMEM GMST signal (method by Kajtar et al. (2019)) at each grid point
        b) remove MMEM SST index (Steinman et al. (2015))

        1. load raw SST data
        2. generate forced signal
            model:  fit to GMST
                linear
                quadratic
            observations:
                single-factor CMIP GMST MMEM
                two-factor CMIP all natural and CMIP anthropogenic (= all forcings - all natural)
        3. regression:
            single time series: forced signal onto SST data -> \beta
            two time series:
        4. use regression coefficient \beta to generate SST signal due to forcing
        5. remove that signal

        run            .. CESM simulation name
        tres           .. time resolution
        detrend_signal .. either GMST (Kajtar et al. (2019))
                          or target region (Steinman et al. (2015))
        time_slice     .. time range selected
        """
        assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'had']
        assert tres in ['yrly', 'monthly']
        assert detrend_signal in ['GMST', 'AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']
        if detrend_signal in ['AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']:
            assert run=='had'
        if run=='had':
            assert time_slice=='full'

        # file name and domain
        fn = f'{path_samoc}/SST/SST_{tres}_{run}.nc'
        if run in ['ctrl', 'rcp']:
            if tres=='yrly':
                domain = 'ocn'
            elif tres=='monthly':
                domain = 'ocn_rect'
        elif run in ['lpd', 'lpi']:  
            domain = 'ocn_low'
        elif run=='had':
            domain = 'ocn_had'

        print('load and subselect data')
        MASK = boolean_mask(domain=domain, mask_nr=0, rounded=True)
        SST = self.select_time_slice(xr.open_dataarray(f'{path_samoc}/SST/SST_{tres}_{run}.nc',\
                                                  decode_times=False).where(MASK),
                                time_slice)
        
        if time_slice!='full':
            first_year, last_year = time_slice
        
        if tres=='monthly':  # deseasonalize
            for t in range(12):
                SST[t::12,:,:] -= SST[t::12,:,:].mean(dim='time')
        SST = SST - SST.mean(dim='time')

        print('calculate forced signal')
        forced_signal = self.forcing_signal(run=run, tres=tres, detrend_signal=detrend_signal, time_slice=time_slice)

        if detrend_signal=='GMST':
            print('Kajtar et al. (2019) scaled MMM GMST detrending method')
            if time_slice=='full':
                fn = f'{path_samoc}/SST/SST_beta_{detrend_signal}_{tres}_{run}.nc'
            else:
                fn = f'{path_samoc}/SST/SST_beta_{detrend_signal}_{tres}_{run}_{first_year}_{last_year}.nc'
                
            try:
                assert 1==0
                assert os.path.exists(fn)
                print('reusing previously calculated beta!')
                print(f'file exists: {fn}')
                beta = xr.open_dataset(fn).slope
            except:
                if run=='ctrl': SST = SST[40:,:,:]
                beta = ADA().lag_linregress(forced_signal, SST)['slope']
                if run=='had':
                    beta = xr.where(abs(beta)<5, beta, np.median(beta))
                ds = xr.merge([forced_signal, beta])
                ds.to_netcdf(fn)
                
            SST_dt = SST - beta * forced_signal
            SST_dt = SST_dt - SST_dt.mean(dim='time')
            print('test')

            # output name
            if run=='had':
                dt = 'sfdt'  # single factor detrending
            elif run in ['ctrl', 'lpd', 'rcp']:
                dt = 'sqdt'  # scaled quadratic detrending
            else:
                dt = 'sldt'  # scaled linear detrending
            
        elif detrend_signal in ['AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']:
            print('Steinman et al. (2015) method')
            # these indices will be detrended afterwards
            SST_dt = SST - forced_signal
            ds = None
            dt = f'{detrend_signal}dt'

        print('writing output')
        if time_slice=='full':
            fn = f'{path_samoc}/SST/SST_{detrend_signal}_{dt}_{tres}_{run}.nc'
        else:
            fn = f'{path_samoc}/SST/SST_{detrend_signal}_{dt}_{tres}_{run}_{first_year}_{last_year}.nc'
        SST_dt.to_netcdf(fn)
        print(f'detrended {run} SST file written out to:\n{fn}')
        
        # additional two factor detrending for had
        if run=='had' and tres=='yrly':
            self.two_factor_detrending(SST)
            
        return


    def forcing_signal(self, run, tres, detrend_signal, time_slice='full'):
        """ GMST forced component
        run            .. dataset
        tres           .. time resolution
        detrend_signal .. 
        time_slice
        """
        print('creating forcing signal')
        assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'had']
        assert tres in ['yrly', 'monthly']
        assert detrend_signal in ['GMST', 'AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']

        # simulations: linear/quadratic fit to GMST signal
        if run in ['ctrl', 'rcp', 'lpd', 'lpi']:
            assert detrend_signal=='GMST'
            if run=='rcp':  # need actual GMST time series
                forced_signal = xr.open_dataset(f'{path_samoc}/GMST/GMST_{tres}_{run}.nc', decode_times=False).GMST
                forced_signal = xr_quadtrend(forced_signal)
            elif run in ['ctrl', 'lpd']:  # use global mean SST as as proxy
                forced_signal = xr.open_dataarray(f'{path_samoc}/SST/GMSST_{tres}_{run}.nc', decode_times=False)
                if run=='ctrl':  # strong adjustment in the first 40 years
                    forced_signal = forced_signal[40:]
                forced_signal = xr_quadtrend(forced_signal)
            else:  # create mock linear trend 
                times = xr.open_dataarray(f'{path_samoc}/SST/SST_{tres}_{run}.nc', decode_times=False).time
                forced_signal = xr.DataArray(np.linspace(0,1,len(times)), coords={'time': np.sort(times.values)}, dims=('time'))
                forced_signal.name = 'GMST'
                forced_signal.attrs = {'Note':'This is the mock linear trend, simply going from 0 to 1'}
#             if tres=='yrly':
#                 times = forced_signal['time'] + 31 # time coordinates shifted by 31 days (SST saved end of January, GMST beginning)
#                 if run=='ctrl':  # for this run, sometimes 31 days, sometimes 15/16 days offset
#                     times = xr.open_dataset(f'{path_samoc}/SST/SST_yrly_ctrl.nc', decode_times=False).time
#                 forced_signal = forced_signal.assign_coords(time=times)

        # observations: CMIP5 multi model ensemble mean of all forcings GMST
        elif run=='had':
            forced_signal = xr.open_dataarray(f'{path_data}/CMIP5/KNMI_CMIP5_{detrend_signal}_{tres}.nc', decode_times=False)
            if tres=='monthly':  # deseasonalize
                for t in range(12):
                    forced_signal[t::12] -= forced_signal[t::12].mean(dim='time')

            if tres=='yrly':# select 1870-2018
                times = (forced_signal['time'].astype(int) - 9)*365
                forced_signal = forced_signal.assign_coords(time=times)  # days since 1861
                forced_signal = forced_signal[9:158]
            elif tres=='monthly':
                # ...
                forced_signal = forced_signal[9*12:158*12-1]
                times = xr.open_dataarray(f'{path_samoc}/SST/SST_monthly_had.nc', decode_times=False).time.values
                forced_signal = forced_signal.assign_coords(time=times)

        forced_signal = self.select_time_slice(forced_signal, time_slice)

        forced_signal -= forced_signal.mean()
        forced_signal.name = 'forcing'
        return forced_signal
    
    
    def SST_pointwise_detrending(self, run, tres='yrly', degree=1, time_slice='full'):
        """ calculates the trends of ocean fields

        input:
        SST    .. xr DataArray
        degree .. degree of polynomial to remove

        output:
        SST_dt .. pointwise detrended SST field

        7 secs for lpd run, 40 seconds
        """
        print('detrending SST pointwise')
        assert degree in [1, 2]
        if run in ['ctrl', 'rcp']:
            MASK = boolean_mask('ocn'     , 0)
        elif run in ['lpi', 'lpd']:
            MASK = boolean_mask('ocn_low' , 0)
        (jm, im) = MASK.shape
        fn = f'{path_samoc}/SST/SST_{tres}_{run}.nc'
        SST = self.select_time_slice(xr.open_dataarray(fn, decode_times=False).where(MASK),
                                     time_slice)
        SST   = SST.where(MASK>0).fillna(-9999)

        Nt = SST.values.shape[0]
        A = SST.values.reshape((Nt, im*jm))

        SST_pf = np.polyfit(SST.time, A, degree)
        
        pf0 = A[0,:].copy()
        pf1 = A[0,:].copy()
        pf0 = SST_pf[0,:]
        pf1 = SST_pf[1,:]
        if degree==1:
#             SST_dt = pf0*SST.time - pf1
            detrend_signal = 'linear'
        elif degree==2:
            pf2 = A[0,:].copy()
            pf2 = SST_pf[2,:]
            A_dt = np.expand_dims(SST.time**2             , 1).dot(np.expand_dims(SST_pf[0,:], 0)) \
                 + np.expand_dims(SST.time                , 1).dot(np.expand_dims(SST_pf[1,:], 0)) \
                 + np.expand_dims(np.ones((len(SST.time))), 1).dot(np.expand_dims(SST_pf[2,:], 0))
            detrend_signal = 'quadratic'
        dt = 'pwdt'  # pointwise detrended

        fn_new = f'{path_samoc}/SST/SST_{detrend_signal}_{dt}_yrly_{run}.nc'
        SST_dt = SST.copy()
        SST_dt.values = (A-A_dt).reshape((Nt,jm,im))
        SST_dt.to_netcdf(fn_new)
        print(f'created {fn_new}')
        return
    
    
    def two_factor_detrending(self, SST):
        print(f'additional two factor detrending for had SST')
        # load CMIP5 multi-model means
        forcing_natural = xr.open_dataarray(f'{path_samoc}/GMST/CMIP5_natural.nc', decode_times=False)
        forcing_anthro  = xr.open_dataarray(f'{path_samoc}/GMST/CMIP5_anthro.nc' , decode_times=False)
        forcing_all     = xr.open_dataarray(f'{path_samoc}/GMST/CMIP5_all.nc'    , decode_times=False)

        for forcing in [forcing_natural, forcing_anthro, forcing_all]:
            forcing.coords['time'] = (forcing.time-9)*365

        forcings = forcing_natural.to_dataframe(name='natural').join(
                     [forcing_anthro.to_dataframe( name='anthro'),
                      forcing_all.to_dataframe(name='all')])

        SST_stacked = SST.stack(z=('latitude', 'longitude'))
        ds_anthro   = SST_stacked[0,:].squeeze().copy()
        ds_natural  = SST_stacked[0,:].squeeze().copy()

        # multiple linear regression
        X = sm.add_constant(forcings[['anthro', 'natural']])
        for i, coordinate in enumerate(SST_stacked.z):
            y = SST_stacked[:, i].values
            model = sm.OLS(y, X).fit()
            ds_anthro[i] = model.params['anthro']
            ds_natural[i] = model.params['natural']

        beta_anthro  = ds_anthro .unstack('z')
        beta_natural = ds_natural.unstack('z')
        
        # output
        ds = xr.merge([{'forcing_anthro': forcing_anthro}, {'beta_anthro': beta_anthro}])
        ds.to_netcdf(f'{path_samoc}/SST/SST_beta_anthro_GMST_yrly_had.nc')
        ds = xr.merge([{'forcing_natural': forcing_natural}, {'beta_natural':beta_natural}])
        ds.to_netcdf(f'{path_samoc}/SST/SST_beta_natural_GMST_yrly_had.nc')

        SST_dt = SST - beta_anthro*forcing_anthro - beta_natural*forcing_natural

        dt = 'tfdt'  # two factor detrending
        fn = f'{path_samoc}/SST/SST_GMST_{dt}_yrly_had.nc'
        SST_dt.to_netcdf(fn)
        print(f'detrended had SST file written out to:\n{fn}')
        return
    
    

    
    
    # AUXILIARY FUNCTIONS
    
    
    def select_time_slice(self, data, time_slice):
        """ if time_slice is not `full`, return subselected data
        time_slice .. either `full` or (start_year, end_year) tuple
        """
        assert time_slice=='full' or type(time_slice)==tuple
        if time_slice!='full':  # use subset in time
            time_coords = tuple(365*x+31 for x in time_slice)
            data = data.sel(time=slice(*time_coords))
        else:  # return original data
            pass
        return data
    
    
    def remove_superfluous_files(self, fn):
        print('removing superfluous files')
        for x in glob.glob(fn):
            os.remove(x) 
        
        
    def shift_had(self, da):
        """ shifts lons to [0,360] to make Pacific contiguous """
        return da.assign_coords(longitude=(da.longitude+360)%360).roll(longitude=180, roll_coords=True)

    def focus_data(self, da):
        """ drops data outside rectangle around Pacific """
        if 't_lat' in da.coords:  # ctrl
            lat, lon = 't_lat', 't_lon'
        elif 'nlat' in da.coords:  # lpd
            lat, lon = 'nlat', 'nlon'
        elif 'latitude' in da.coords:  # had
            lat, lon = 'latitude', 'longitude'
        else:  raise ValueError('xr DataArray does not have the right lat/lon coords.')
        da = da.dropna(dim=lat, how='all')
        da = da.dropna(dim=lon, how='all')
        return da


#     def determine_years_from_slice(run, tres, time_slice):
#         assert time_slice is not 'full'
#         if tres=='yrly':
#             first_year, last_year = int(time_slice[0]/365), int(time_slice[1]/365)
#         elif tres=='monthly':
#             if run in ['ctrl', 'rcp']:
#                 first_year, last_year = int(time_slice[0]/12), int(time_slice[1]/12)
#                 if run=='ctrl':
#                     first_year += 100
#                     last_year  += 100
#                 elif run=='rcp':
#                     first_year += 2000
#                     last_year  += 2000
#             elif run in ['lpd', 'lpi']:
#                 first_year, last_year = int(time_slice[0]/365), int(time_slice[1]/365)
#         return first_year, last_year