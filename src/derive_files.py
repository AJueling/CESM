import os
import xarray as xr

from paths import path_samoc
from timeseries import IterateOutputCESM
from xr_DataArrays import depth_lat_lon_names, xr_DZ, xr_DXU
from xr_regression import xr_autocorrelation_2D, lag_linregress_3D


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
        print(f'self.domain = {self.domain} for {self.run} run, some functions require setting different domain')
    
    def multifile_name(self)
    
    
    def yrly_avg_nc(self, fields, test=False):
        """ creates yearly average file from monthly
        input:    fields .. list of field names
        output:   [writes netCDF file]
        (takes approx. 2 min for high res ocean data for two 3D fields)
        (takes approx. 4 sec for lower res atm data for one 2D and one 3D field)
        """
        assert self.domain in ['ocn', 'ocn_rect', 'atm', 'ice']
        assert self.run in ['ctrl', 'rcp', 'lpd' ,'lpi']

        print(f'yearly averaging of {self.run} {self.domain}')
        for field in fields:  print(f'   {field}')

        name = ''
        n_fields = len(fields)
        for i, field in enumerate(fields):
            name += field
            if i<n_fields-1:
                name += '_'

        ffield = fields[0]

        first_year = IterateOutputCESM(domain=self.domain, run=self.run, tavg='monthly').year

        for y, m, s in IterateOutputCESM(domain=self.domain, run=self.run, tavg='monthly'):
            if m==1:
                new_filename = CESM_filename(domain=self.domain, run=self.run, y=y, m=0, name=name)
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
    
    def make_detrended_SST_file(self, )
    
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
            for i, (y,m,s) in enumerate(IterateOutputCESM(domain=self.domain, run=self.run, tavg='yrly', name='UVEL_VVEL')):
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
