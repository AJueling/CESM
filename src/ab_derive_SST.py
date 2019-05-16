import os
import numpy as np
import xarray as xr

from paths import CESM_filename, path_samoc
from timeseries import IterateOutputCESM
from xr_DataArrays import depth_lat_lon_names, xr_DZ, xr_DXU


class DeriveSST(object):
    """ generate fields """
    def __init__(self):
        return
    
    @staticmethod
    def remove_superfluous_files(fn):
        for x in glob.glob(fn):
            os.remove(x) 
    
    @staticmethod
    def generate_yrly_SST_files(run):
        """ generate the SST data files from TEMP_PD yearly averaged files """
        # ca. 4:30 min for ctrl/rcp, 1:25 for lpi
        # stacking files into one xr DataArray object
        for i, (y,m,s) in enumerate(IterateOutputCESM('ocn', run, 'yrly', name='TEMP_PD')):
            print(y)
            da = xr.open_dataset(s, decode_times=False).TEMP[0,:,:]
            da = da.drop(['z_t', 'ULONG', 'ULAT'])
            da['TLAT' ] = da['TLAT' ].round(decimals=2)
            da['TLONG'] = da['TLONG'].round(decimals=2)
            del da.encoding["contiguous"]
            ds = t2ds(da=da, name='SST', t=int(round(da.time.item())))
            ds.to_netcdf(path=f'{path_samoc}/SST/SST_yrly_{run}_{y}.nc', mode='w')

        combined = xr.open_mfdataset(f'{path_samoc}/SST/SST_yrly_{run}_*.nc',
                                     concat_dim='time', autoclose=True, coords='minimal')
        combined.to_netcdf(f'{path_samoc}/SST/SST_yrly_{run}.nc')

        GenerateSSTFields.remove_superfluous_files(f'{path_samoc}/SST/SST_yrly_{run}_*.nc')
            
            
    @staticmethod
    def generate_monthly_SST_files(run):
        """ concatonate monthly files, ocn_rect for high res runs"""
        # 8 mins for 200 years of ctrl
        if run in ['ctrl', 'rcp']:  domain = 'ocn_rect'
        elif run in ['lpd', 'lpi']:  domain = 'ocn_low'
            
        for y,m,s in IterateOutputCESM(domain=domain, tavg='monthly', run=run):
            if run in ['ctrl', 'rcp']:
                xa = xr.open_dataset(s, decode_times=False).TEMP[0,:,:]
            if run in ['lpd', 'lpi']:
                xa = xr.open_dataset(s, decode_times=False).TEMP[0,0,:,:]
            if m==1:
                print(y)
                xa_out = xa.copy()    
            else:
                xa_out = xr.concat([xa_out, xa], dim='time')
            if m==12:
                xa_out.to_netcdf(f'{path_samoc}/SST/SST_monthly_{run}_{y}.nc')
                        
        combined = xr.open_mfdataset(f'{path_samoc}/SST/SST_monthly_{run}_*.nc',
                                     concat_dim='time', decode_times=False)
        combined.to_netcdf(f'{path_samoc}/SST/SST_monthly_{run}.nc')
        combined.close()

        GenerateSSTFields.remove_superfluous_files(f'{path_samoc}/SST/SST_monthly_{run}_*.nc')
            
            
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
                    xa_out.to_netcdf(f'{path_samoc}/SST/SST_monthly_{r}_{run}_{y}.nc')
                            
            combined = xr.open_mfdataset(f'{path_samoc}/SST/SST_monthly_{r}_{run}_*.nc',
                                         concat_dim='time', decode_times=False)
            combined.to_netcdf(f'{path_samoc}/SST/SST_monthly_{r}_{run}.nc')
            combined.close()
            
            GenerateSSTFields.remove_superfluous_files(f'{path_samoc}/SST/SST_monthly_{r}_{run}_*.nc')
            # # remove yearly files
            
    @staticmethod
    def generate_monthly_mock_linear_GMST_files(run):
        """ generates xr DataArray with same time coordinates as monthly SST fields
        which contains the linear fit to the mean SST as a stand in for the missing
        monthly GMST
        """
        assert run in ['ctrl', 'lpd']
        if run=='ctrl':  dims = ['t_lon', 't_lat']
        elif run=='lpd':  dims = ['nlon', 'nlat']
        da = xr.open_dataarray(f'{path_samoc}/SST/SST_monthly_{run}.nc', decode_times=False)
        da_new = xr_lintrend(da.mean(dim=dims, skipna=True, keep_attrs=True))
        da_new.name = 'GMST'
        da_new.attrs = {'Note':'This is the linear trend of the SST evolution, not GMST'}
        da_new.to_dataset().to_netcdf(f'{path_samoc}/GMST/GMST_monthly_{run}.nc')
        
        
    @staticmethod
    def SST_remove_forced_signal(run, tres='yrly', detrend_signal='GMST', time_slice='full'):
        """ removed the scaled, forced GMST signal (method by Kajtar et al. (2019))

        1. load raw SST data
        2. generate forced signal (either quadtrend or CMIP MMEM)
        3. regress forced signal onto SST data -> \beta
        4. use regression coefficient \beta to generate SST signal due to forcing
        5. remove that signal

        run            ..
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

        first_year, last_year = self.determine_years_from_slice(run, tres, time_slice)

    #     sys.exit("Error message")
        # 1. load data
        MASK = boolean_mask(domain=domain, mask_nr=0, rounded=True)
        SST = xr.open_dataarray(f'{path_samoc}/SST/SST_{tres}_{run}.nc', decode_times=False).where(MASK)
        if time_slice is not 'full':
            assert type(time_slice)==tuple
            SST = SST.sel(time=slice(*time_slice))

        if tres=='monthly':  # deseasonalize
            for t in range(12):
                SST[t::12,:,:] -= SST[t::12,:,:].mean(dim='time')

        SST = SST - SST.mean(dim='time')

        # 2/3/4. calculate forced signal
        forced_signal = self.forcing_signal(run=run, tres=tres, detrend_signal=detrend_signal, time_slice=time_slice)

        if detrend_signal=='GMST':
            beta = lag_linregress_3D(forced_signal, SST)['slope']
            if run=='had':
                beta = xr.where(abs(beta)<5, beta, np.median(beta))
            ds = xr.merge([forced_signal, beta])
            ds.to_netcdf(f'{path_samoc}/SST/SST_beta_{detrend_signal}_{tres}_{run}.nc')
            forced_map = beta * forced_signal

            # 5.
            SST_dt = SST - forced_map
            SST_dt -= SST_dt.mean(dim='time')

        elif detrend_signal in ['AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']:
            # these indices will be detrended afterwards
            SST_dt = SST - forced_signal
            ds = None

        if time_slice=='full':
            SST_dt.to_netcdf(f'{path_samoc}/SST/SST_{detrend_signal}_dt_{tres}_{run}.nc')
        else:
            SST_dt.to_netcdf(f'{path_samoc}/SST/SST_{detrend_signal}_dt_{tres}_{first_year}_{last_year}_{run}.nc')

        return SST_dt, ds


    def forcing_signal(run, tres, detrend_signal, time_slice='full'):
        """ GMST forced component
        run            .. dataset
        tres           .. time resolution
        detrend_signal .. 
        time_slice
        """
        assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'had']
        assert tres in ['yrly', 'monthly']
        assert detrend_signal in ['GMST', 'AMO', 'SOM', 'TPI1', 'TPI2', 'TPI3']

        if run in ['ctrl', 'rcp', 'lpd', 'lpi']:
            assert detrend_signal=='GMST'
            forced_signal = xr.open_dataset(f'{path_samoc}/GMST/GMST_{tres}_{run}.nc', decode_times=False).GMST
            if run=='rcp':
                forced_signal = xr_quadtrend(forced_signal)
            else:
                forced_signal = xr_lintrend(forced_signal)
            if tres=='yrly':
                times = forced_signal['time'] + 31 # time coordinates shifted by 31 days (SST saved end of January, GMST beginning)
                if run=='ctrl':  # for this run, sometimes 31 days, sometimes 15/16 days offset
                    times = xr.open_dataset(f'{path_samoc}/SST/SST_yrly_ctrl.nc', decode_times=False).time
                forced_signal = forced_signal.assign_coords(time=times)
    #         elif tres=='monthly':

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

        if time_slice is not 'full':
            forced_signal = forced_signal.sel(time=slice(*time_slice))

        forced_signal -= forced_signal.mean()
        forced_signal.name = 'forcing'

        return forced_signal
        
        
    def determine_years_from_slice(run, tres, time_slice):
        assert time_slice is not 'full'
        if tres=='yrly':
            first_year, last_year = int(time_slice[0]/365), int(time_slice[1]/365)
        elif tres=='monthly':
            if run in ['ctrl', 'rcp']:
                first_year, last_year = int(time_slice[0]/12), int(time_slice[1]/12)
                if run=='ctrl':
                    first_year += 100
                    last_year  += 100
                elif run=='rcp':
                    first_year += 2000
                    last_year  += 2000
            elif run in ['lpd', 'lpi']:
                first_year, last_year = int(time_slice[0]/365), int(time_slice[1]/365)
        return first_year, last_year