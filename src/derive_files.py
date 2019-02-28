import os
import xarray as xr

from paths import path_samoc
from timeseries import IterateOutputCESM
from xr_DataArrays import depth_lat_lon_names
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
    
    def make_GMST_file(self):
    
    
    def make_OHC_file(self):
        return
    
    def make_MOC_file(self):
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
