import os
import numpy as np
import xarray as xr
import scipy.signal as signal

from paths import path_ocn_ctrl, path_ocn_rcp
from paths import path_atm_ctrl, path_atm_rcp
from paths import path_yrly_ctrl, path_yrly_rcp
from paths import rcpstr, spinup, CESM_filename
from paths import path_samoc, path_results
from constants import abs_zero
from xr_integrate import xr_surf_mean, xr_zonal_mean
from xr_DataArrays import xr_AREA


class IterateOutputCESM:
    """ iterator over all CESM ctrl/rcp filenames
    automatically detects the last file
    
    example:
    >for year, month, filename in IterateOutputCESM('ocn', ctrl', 'monthly'):
    >    print(year, month, filename)
    """
    
    def __init__(self, domain, run, tavg, name=None):
        assert domain in ['ocn', 'ocn_rect', 'ocn_low', 'atm', 'ice']
        assert run in ['ctrl', 'rcp', 'lpd', 'lpi']
        assert tavg in ['monthly', 'yrly']
        
        self.domain = domain
        self.run    = run
        self.tavg   = tavg
        self.stop   = False
        self.name   = name
        
        if tavg=='monthly':  self.month = 1
        elif tavg=='yrly':   self.month = 0
        
        if run=='ctrl':   self.year  =  100
        elif run=='rcp':  self.year  = 2000
        elif run=='lpd':  self.year  =  154
        elif run=='lpi':  
            if domain in ['ocn', 'ocn_low']:
                self.year  = 1600
            if domain=='atm':
                self.year  = 2876
            
    def file(self):
        if self.tavg=='monthly':
            filename = CESM_filename(self.domain, self.run, self.year, self.month)     
        elif self.tavg=='yrly':
            if self.domain=='ocn':
                if self.name==None:
                    raise ValueError('must provide (variables part of) name for yrly file')
                else:
                    filename = CESM_filename(self.domain, self.run, self.year, self.month, self.name)
            elif self.domain=='atm':
                if self.run in ['ctrl', 'rcp', 'lpi']:
                    filename = CESM_filename(self.domain, self.run, self.year, self.month, self.name)
                elif self.run in ['lpd']:  # yrly files are written out already
                    if self.name!=None:  print("name is ignored, as yearly files existed already")
                    filename = CESM_filename(domain=self.domain, run=self.run, y=self.year, m=self.month, name=None)
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
            if self.tavg=='monthly':
                self.month += 1
                length +=1
                if self.month==13:
                    self.month = 1
                    self.year +=1
            elif self.tavg=='yrly':
                self.year +=1
                length +=1
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
            


def ncfile_list(domain, run, tavg, name=None):
    list_files = []
    for (y,m,s) in IterateOutputCESM(domain=domain, run=run, tavg=tavg, name=name):
        assert os.path.exists(s)
        list_files.append(s)
    return list_files


###############################################################################
# time averaging:  monthly -> yearly
###############################################################################


def yrly_avg_nc(domain, run, fields, test=False):
    """ creates yearly average file
    
    input:
    domain .. (str) 'ocn' or 'atm'
    run    .. (str) 'ctrl', 'rcp', 'lpd' ,'lpi'
    fields .. list of field names
    
    output:
    [writes netCDF file]
    
    (takes approx. 2 min for high res ocean data for two 3D fields)
    (takes approx. 4 sec for lower res atm data for one 2D and one 3D field)
    """
    assert domain in ['ocn', 'ocn_rect', 'atm', 'ice']
    assert run in ['ctrl', 'rcp', 'lpd' ,'lpi']
    
    print(f'yearly averaging of {run} {domain}')
    for field in fields:  print(f'   {field}')
    
    name = ''
    n_fields = len(fields)
    for i, field in enumerate(fields):
        name += field
        if i<n_fields-1:
            name += '_'
    
    ffield = fields[0]
    
    first_year = IterateOutputCESM(domain=domain, run=run, tavg='monthly').year
    
    for y, m, s in IterateOutputCESM(domain=domain, run=run, tavg='monthly'):
        if m==1:
            new_filename = CESM_filename(domain=domain, run=run, y=y, m=0, name=name)
        if os.path.exists(new_filename):  continue
        
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


###############################################################################
# time filter
###############################################################################


def lowpass(ts, period):
    """ lowpass filter
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.butter.html
    ts      .. time series array with axis 0 as time axis
    perdiod .. cutoff period
    """
    N  = 2         # Filter order
    Wn = 1/period  # Cutoff frequency
    B, A = signal.butter(N, Wn, output='ba')
    filtered = signal.filtfilt(B, A, ts, axis=(0), padlen=N-1, padtype='constant')
    
    if type(ts)==xr.core.dataarray.DataArray:
        ts_new = ts.copy()
        ts_new.values = filtered
    else:
        ts_new = filtered
    
    return ts_new


def chebychev(ts, period):
    N, rp = 6, 1  # N=6 was used in Henley et al. (2015), rp is low to minimize ripples
    Wn = 1/period
    B, A = signal.cheby1(N, rp, Wn)
    filtered = signal.filtfilt(B, A, ts, axis=(0), padlen=N-1, padtype='constant')
    
    if type(ts)==xr.core.dataarray.DataArray:
        ts_new = ts.copy()
        ts_new.values = filtered
    else:
        ts_new = filtered
    
    return ts_new


def notch(ts, period):
    """ single frequency filter 
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.iirnotch.html
    period .. 1/f to filter out
    """
    w0 = 1/period
    Q  = .2# found to work fairly well for deseasonalizing
    B, A = signal.iirnotch(w0, Q)
    filtered = signal.filtfilt(B, A, ts, axis=(0), padlen=12)
    
    if type(ts)==xr.core.dataarray.DataArray:
        ts_new = ts.copy()
        ts_new.values = filtered
    else:
        ts_new = filtered
    
    return ts_new