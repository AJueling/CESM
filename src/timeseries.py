import os
import numpy as np
import xarray as xr

from paths import path_ocn_ctrl, path_ocn_rcp
from paths import path_atm_ctrl, path_atm_rcp
from paths import path_yrly_ctrl, path_yrly_rcp
from paths import rcpstr, spinup, CESM_filename
from paths import path_samoc
from constants import abs_zero
from xr_integrate import xr_surf_mean, xr_zonal_mean
from xr_DataArrays import generate_xr_AREA

class IterateOutputCESM:
    """ iterator over all CESM ctrl/rcp filenames
    automatically detects the last file
    
    example:
    >for filename in IterateOutputCESM('ctrl'):
    >    print(filename)
    """
    
    def __init__(self, domain, run, tavg):
        assert domain=='ocn' or domain=='atm'
        assert run=='rcp' or run=='ctrl'
        
        self.domain = domain
        self.run    = run
        self.tavg   = tavg
        self.stop   = False
        
        if tavg=='monthly':  self.month  = 1
        elif tavg=='yrly':   self.month  = 0
        
        if run=='rcp':     self.year  = 2000
        elif run=='ctrl':  self.year  = 200
            
    def file(self):
        if self.tavg=='monthly':
            filename = CESM_filename(self.domain, self.run, self.year, self.month)     
        elif self.tavg=='yrly':
            if self.domain=='atm':
                filename = CESM_filename(self.domain, self.run, self.year, self.month, name='T_T850')  
        if os.path.exists(filename)==False:
            self.stop = True
        return filename
    
    def __iter__(self):
        return self
    
    def __len__(self):
        length = 0
        year   = self.year
        month  = self.month
            
        while os.path.exists(self.file()):
            self.month += 1
            if self.tavg=='monthly':  length +=1
            if self.month==13:
                self.month = 1
                self.year +=1
                if self.tavg=='yrly':  length +=1
        return length

    def __next__(self):
        new_file = self.file()
        y = self.year
        m = self.month
        if self.stop:
            raise StopIteration
        else:
            if self.tavg=='monthly':
                if self.month==12:
                    self.month = 1
                    self.year +=1
                else:
                    self.month += 1
                return y, m, new_file
            
            elif self.tavg=='yrly':
                self.year += 1
                return y, 0, new_file



def yrly_avg_nc(domain, run, fields, test=False):
    """ creates yearly average file
    
    input:
    domain .. (str) 'ocn' or 'atm'
    run    .. (str) 'ctrl' or 'rcp'
    fields .. list of field names
    
    output:
    [writes netCDF file]
    
    takes approx. 2 min for high res ocean data for two 3D fields
    takes approx. 4 sec for lower res atm data for one 2D and one 3D field 
    """
    assert domain=='ocn' or domain=='atm'
    assert run=='rcp' or run=='ctrl'
    
    name = ''
    n_fields = len(fields)
    for i, field in enumerate(fields):
        name += field
        if i<n_fields-1:
            name += '_'
    
    ffield = fields[0]
    i = -1
    
    for y, m, s in IterateOutputCESM(domain=domain, run=run, tavg='monthly'):
        
        print(m)
        ds = xr.open_dataset(s, decode_times=False)
        
        if m==1:  # create new xr Dataset
            i += 1
            dim = len(np.shape(ds[ffield]))
            if dim==3:  # 2D field
                ds_out = (ds[ffield][0,:,:]/12).to_dataset()
            elif dim==4:  # 3D
                ds_out = (ds[ffield][0,:,:,:]/12).to_dataset()
                
            for field in fields[1:]:  # add rest of fields
                dim = len(np.shape(ds[field]))
                if dim==3:
                    ds_out[field] = ds[field][0,:,:]/12
                elif dim==4:
                    ds_out[field] = ds[field][0,:,:,:]/12
        else:
            for field in fields:
                dim = len(np.shape(ds[field]))
                if   dim==3:  ds_out[field][:,:]   += ds[field][0,:,:]/12
                elif dim==4:  ds_out[field][:,:,:] += ds[field][0,:,:,:]/12

        if m==12 and test==False:  # write to new file
            print(y, CESM_filename(domain, run, y, 0))
            ds_out.to_netcdf(path=CESM_filename(domain, run, y, 0, name=name),
                             mode='w')
            
        if test==True and y==IterateOutputCESM(domain, run).year+1:
            break
    
    return



def GMST_timeseries(run):
    """ builds a timesries of the GMST and saves it to a netCDF
    
    input:
    run .. (str) ctrl or cp
    
    output:
    ds_new  .. xr Dataset containing GMST and T_zonal
    """
    
    dx = 1
    nlats = int(180/dx)
    ny = len(IterateOutputCESM(domain='atm', run=run, tavg='yrly'))
    years = np.arange(ny) + IterateOutputCESM(domain='atm', run=run, tavg='yrly').year
    lats = np.arange(-90+dx/2, 90, dx)
    assert len(lats)==nlats
    
    AREA = generate_xr_AREA('atm')
    
    for i, (y, m, file) in enumerate(IterateOutputCESM(domain='atm', run=run, tavg='yrly')):
        ds = xr.open_dataset(file, decode_times=False)
        if i==0:
            ds_new = xr.Dataset()
            
            da1 = xr.DataArray(data=np.zeros((ny)),
                               coords={'year': years},
                               dims=('year'))
            da2 = xr.DataArray(data=np.zeros((ny, nlats)),
                               coords={'year': years, 'lat_bins': lats},
                               dims=('year', 'lat_bins'))
            ds_new['GMST']    = da1
            ds_new['T_zonal'] = da2
            
        ds_new['GMST'][i]      = xr_surf_mean(ds['T'][-1,:,:], AREA=AREA) + abs_zero
        ds_new['T_zonal'][i,:] = xr_zonal_mean(ds['T'][-1,:,:], AREA=AREA, dx=1, lat_name='lat') + abs_zero
        
    ds_new.to_netcdf(path=f'{path_samoc}/{run}/GMST.nc', mode='w')
    
    return ds_new