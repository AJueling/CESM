import os
import numpy as np
import xarray as xr

from paths import CESM_filename, path_samoc
from timeseries import IterateOutputCESM
from xr_DataArrays import depth_lat_lon_names, xr_DZ, xr_DXU
# from xr_regression import xr_autocorrelation_2D, lag_linregress_3D


class MakeDerivedFiles(object):
    """ functions to generate netcdf files derived from CESM output / obs. """
    
    def __init__(self, run):
        assert run in ['ctrl', 'rcp', 'lpd', 'lpi', 'had']
        self.run = run
        if run in ['ctrl', 'rcp']:
            self.domain = 'ocn'
        elif run in ['lpd', 'lpi']:
            self.domain = 'ocn_low'
        elif run=='had':
            self.domain = 'ocn_had'
        print(f'self.domain = {self.domain} for {self.run} run, some functions may require setting different domain')
    
    
    def multifile_name(self):
        return
    
    
    def yrly_avg_nc(self, domain, fields, test=False):
        """ creates yearly average file from monthly
        input:    fields .. list of field names
        output:   [writes netCDF file]
        (takes approx. 2 min for high res ocean data for two 3D fields)
        (takes approx. 4 sec for lower res atm data for one 2D and one 3D field)
        """
        assert domain in ['ocn', 'ocn_rect', 'atm', 'ice']
        assert self.run in ['ctrl', 'rcp', 'lpd' ,'lpi']

        print(f'yearly averaging of {self.run} {domain}')
        for field in fields:  print(f'   {field}')

        name = ''
        n_fields = len(fields)
        for i, field in enumerate(fields):
            name += field
            if i<n_fields-1:
                name += '_'
                
        ffield = fields[0]
                
        if self.run=='lpd' and domain=='atm':
            # there are no monthly files for the lpd run
            # but for consistency, I create the same files as for the runs with monthly output only
            for y, m, s in IterateOutputCESM(domain=domain, run=self.run, tavg='yrly'):
                new_filename = CESM_filename(domain=domain, run=self.run, y=y, m=0, name=name)
                ds = xr.open_dataset(s, decode_times=False)
                dim = len(np.shape(ds[ffield]))
                if   dim==3:  ds_out = (ds[ffield][0,:,:]).to_dataset()
                elif dim==4:  ds_out = (ds[ffield][0,:,:,:]).to_dataset()
                if len(fields)>1:
                    for field in fields[1:]:
                        dim = len(np.shape(ds[field]))
                        if   dim==3:  ds_out[field] = ds[field][0,:,:]
                        elif dim==4:  ds_out[field] = ds[field][0,:,:,:]
                print(y, new_filename)
                ds_out.to_netcdf(path=new_filename, mode='w')
            return
        
        else:                
            first_year = IterateOutputCESM(domain=domain, run=self.run, tavg='monthly').year

            for y, m, s in IterateOutputCESM(domain=domain, run=self.run, tavg='monthly'):
                if m==1:
                    new_filename = CESM_filename(domain=domain, run=self.run, y=y, m=0, name=name)
                if os.path.exists(new_filename):
                    continue
                ds = xr.open_dataset(s, decode_times=False)

                if m==1:  # create new xr Dataset
                    dim = len(np.shape(ds[ffield]))
                    if domain in ['atm', 'ocn', 'ice']:
                        if dim==3:  # 2D field
                            ds_out = (ds[ffield][0,:,:]/12).to_dataset()
                        elif dim==4:  # 3D
                            ds_out = (ds[ffield][0,:,:,:]/12).to_dataset()
                    elif domain=='ocn_rect':
                        if dim==2:  # 2D
                            ds_out = (ds[ffield][:,:]/12).to_dataset()
                        elif dim==3:  # 3D
                            ds_out = (ds[ffield][:,:,:]/12).to_dataset()
                    for field in fields[1:]:  # add rest of fields
                        dim = len(np.shape(ds[field]))
                        if domain in ['atm', 'ocn', 'ice']:
                            if dim==3:
                                ds_out[field] = ds[field][0,:,:]/12
                            elif dim==4:
                                ds_out[field] = ds[field][0,:,:,:]/12
                        elif domain=='ocn_rect':
                            if dim==2:
                                ds_out[field] = ds[field][:,:]/12
                            elif dim==3:
                                ds_out[field] = ds[field][:,:,:]/12

                else:  # add subsequent monthly values
                    for field in fields:
                        dim = len(np.shape(ds[field]))
                        if domain in ['atm', 'ocn', 'ice']:
                            if   dim==3:  ds_out[field][:,:]   += ds[field][0,:,:]/12
                            elif dim==4:  ds_out[field][:,:,:] += ds[field][0,:,:,:]/12
                        elif domain=='ocn_rect':
                            if   dim==2:  ds_out[field][:,:]   += ds[field][:,:]/12
                            elif dim==3:  ds_out[field][:,:,:] += ds[field][:,:,:]/12

                if m==12:  # write to new file
                    print(y, new_filename)
                    ds_out.to_netcdf(path=new_filename, mode='w')

                if test==True and y==first_year+2:  break
            print('done!')
            return
    
    
    def make_SST_yrly_data_file(self):
        """ generates annual SST .nc file (t,lat,lon)
        ca. 4:30 min for ctrl/rcp, 1:25 for lpi
        """
        for i, (y,m,s) in enumerate(IterateOutputCESM('ocn', self.run, 'yrly', name='TEMP_PD')):
            print(y)
            da = xr.open_dataset(s, decode_times=False).TEMP[0,:,:]
            da = da.drop(['z_t', 'ULONG', 'ULAT'])
            da['TLAT' ] = da['TLAT' ].round(decimals=2)
            da['TLONG'] = da['TLONG'].round(decimals=2)
            del da.encoding["contiguous"]
            ds = t2ds(da=da, name='SST', t=int(round(da.time.item())))
            ds.to_netcdf(path=f'{path_samoc}/SST/SST_yrly_{self.run}_{y}.nc', mode='w')

        combined = xr.open_mfdataset(f'{path_samoc}/SST/SST_yrly_{self.run}_*.nc',
                                     concat_dim='time',
                                     autoclose=True,
                                     coords='minimal')
        combined.to_netcdf(f'{path_samoc}/SST/SST_yrly_{self.run}.nc')
        # remove extra netCDF files
        return
    
    
    def make_detrended_SST_file(self):
        return
    
    
    def make_GMST_with_trends_file(self):
        """ builds a timesries of the GMST and saves it to a netCDF
    
        input:
        run    .. (str) ctrl or cp

        output:
        ds_new .. xr Dataset containing GMST and T_zonal
        """
        domain = 'atm'
        tavg   = 'yrly'
        name   = 'T_T850_U_V'

        if run in ['ctrl', 'rcp']:   AREA = xr_AREA('atm')
        elif run=='lpi':             AREA = xr_AREA('atm_f19')
        elif run=='lpd':             AREA = xr_AREA('atm_f09')

        AREA_lat   = AREA.sum(dim='lon')
        AREA_total = AREA.sum(dim=('lat','lon'))


        if run in ['lpd']:  name = None

        ny   = len(IterateOutputCESM(domain=domain, run=run, tavg=tavg, name=name))
        first_yr = IterateOutputCESM(domain=domain, run=run, tavg=tavg, name=name).year
        iterator = IterateOutputCESM(domain=domain, run=run, tavg=tavg, name=name)
        years    = (np.arange(ny) + first_yr)*365  # this is consistent with CESM output

        for i, (y, m, file) in enumerate(iterator):
            print(y)
            assert os.path.exists(file)
            if run in ['ctrl', 'rcp', 'lpi']:
                da = xr.open_dataset(file, decode_times=False)['T'][-1,:,:]
            elif run in ['lpd']:
                da = xr.open_dataset(file, decode_times=False)['T'][0,-1,:,:]

            if i==0:  # create new xr Dataset
                lats = da.lat.values
                ds_new = xr.Dataset()
                ds_new['GMST']    = xr.DataArray(data=np.zeros((ny)),
                                                 coords={'time': years},
                                                 dims=('time'))
                ds_new['T_zonal'] = xr.DataArray(data=np.zeros((ny, len(lats))), 
                                                 coords={'time': years, 'lat': lats},
                                                 dims=('time', 'lat'))

            ds_new['GMST'][i]      = (da*AREA).sum(dim=('lat','lon'))/AREA_total
            ds_new['T_zonal'][i,:] = (da*AREA).sum(dim='lon')/AREA_lat

        # [K] to [degC]
        for field in ['GMST', 'T_zonal']:  
            ds_new[field] = ds_new[field] + abs_zero

        # rolling linear trends [degC/yr]
        ds_new = rolling_lin_trends(ds=ds_new, ny=ny, years=years)

        # fits
        lfit = np.polyfit(np.arange(ny), ds_new.GMST, 1)
        qfit = np.polyfit(np.arange(ny), ds_new.GMST, 2)

        ds_new[f'lin_fit']  = xr.DataArray(data=np.empty((len(ds_new['GMST']))),
                                           coords={'time': years},
                                           dims=('time'),
                                           attrs={'lin_fit_params':lfit})
        ds_new[f'quad_fit'] = xr.DataArray(data=np.empty((len(ds_new['GMST']))),
                                           coords={'time': years},
                                           dims=('time'),
                                           attrs={'quad_fit_params':qfit})

        for t in range(ny):
            ds_new[f'lin_fit'][t]  =                lfit[0]*t + lfit[1]
            ds_new[f'quad_fit'][t] = qfit[0]*t**2 + qfit[1]*t + qfit[2]

        ds_new.to_netcdf(path=f'{path_results}/GMST/GMST_{run}.nc', mode='w')

        return ds_new
        
        
    def make_AHC_file(self):
        """ calculate total atmospheric heat content [Joule] """
        def atm_heat_content(ds, dp, AREA):
            assert 'T' in ds
            assert 'lat' in ds.coords
            assert 'lat' in AREA.coords
            assert 'lon' in ds.coords
            assert 'lon' in AREA.coords
            assert 'lev' in ds.coords
            assert 'lev' in dp.coords
            return (AREA*dp*ds['T']*cp_air).sum()
    
    
    def make_OHC_file(self):
        return
    
    
    def make_MOC_file(self):
        # # ca 20 sec per file
        # # 50 min for ctrl
        # # 26 min for rcp
        if selfrun in ['ctrl', 'rcp']:
            DXU  = xr_DXU(self.domain)            # [m]
            DZU  = xr_DZ(self.domain, grid='U')   # [m]
            # MASK = Atlantic_mask('ocn')   # Atlantic
            # MASK = boolean_mask('ocn', 2) # Pacific
            MASK = boolean_mask(self.domain, 0)   # Global OCean
            for i, (y,m,s) in enumerate(IterateOutputCESM(domain=self.domain,
                                                          run=self.run,
                                                          tavg='yrly',
                                                          name='UVEL_VVEL')):
                # ca. 20 sec per year
                print(i, y, s)
                ds = xr.open_dataset(s, decode_times=False)
                MOC = calculate_MOC(ds=ds, DXU=DXU, DZU=DZU, MASK=MASK)
                if i==0:
                    MOC_out = MOC.copy()
                else:
                    MOC_out = xr.concat([MOC_out, MOC], dim='time')
            #     if y==202: break
            MOC_out.to_netcdf(f'{path_results}/MOC/GMOC_{self.run}.nc')
        
        elif self.run in ['lpd', 'lpi']:
            ds = xr.open_mfdataset()
        
        return
    
    def make_CICE_file(self):
        return
    
    def make_freshwater_35S_file(self):
        return
    
    def make_SST_index_file(self):
        return
    
#     def load_SST_data(self, detrend=False):
#         """ loads raw or detrended SST fields """
#         if detrend==False:
#             self.SST = xr.open_dataarray(f'{path_samoc}/SST/SST_yrly_{self.run}.nc', decode_times=False)
#         if detrend==True:
#             self.SST = xr.open_dataarray(f'{path_samoc}/SST/SST_GMST_dt_yrly_{self.run}.nc', decode_times=False)
#         return



class GenerateSSTFields(object):
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