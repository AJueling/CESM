import os
import numpy as np
import xarray as xr

from paths import path_ocn_ctrl, path_ocn_rcp
from paths import path_atm_ctrl, path_atm_rcp
from paths import path_yrly_ctrl, path_yrly_rcp
from paths import rcpstr, spinup, CESM_filename

class IterateOutputCESM:
    """ iterator over all CESM ctrl/rcp filenames
    automatically detects the last file
    
    example:
    >for filename in IterateOutputCESM('ctrl'):
    >    print(filename)
    """
    
    def __init__(self, domain, run):
        assert domain=='ocn' or domain=='atm'
        assert run=='rcp' or run=='ctrl'
        
        self.domain = domain
        self.run    = run
        self.stop   = False
        self.month  = 1
        
        if run=='rcp':
            self.year  = 2000
            if domain=='ocn':
                self.path  = path_ocn_rcp
            elif domain=='atm':
                self.path  = path_atm_rcp
        elif run=='ctrl':
            self.year  = 200
            if domain=='ocn':
                self.path  = path_ocn_ctrl
            elif domain=='atm':
                self.path  = path_atm_ctrl
            
    def file(self):
        filename = CESM_filename(self.domain, self.run, self.year, self.month)     
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
            if self.month==13:
                self.month = 1
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
            if self.month==12:
                self.month = 1
                self.year +=1
            else:
                self.month += 1
            return y, m, new_file



def yrly_avg_nc(domain, run, fields):
    """ creates yearly average file
    
    input:
    domain .. (str) 'ocn' or 'atm'
    run    .. (str) 'ctrl' or 'rcp'
    fields .. list of field names
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
    
    for y, m, s in IterateOutputCESM(domain, run):
        
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
            
        for field in fields:
            dim = len(np.shape(ds[field]))
            if   dim==3:  ds_out[field][:,:]   += ds[field][0,:,:]/12
            elif dim==4:  ds_out[field][:,:,:] += ds[field][0,:,:,:]/12

        if m==12:  # write to new file
            print(y, CESM_filename(domain, run, y, 0))
            ds_out.to_netcdf(path=CESM_filename(domain, run, y, 0, name=name),
                             mode='w')
            
        if y==IterateOutputCESM(domain, run).year+1:  # for testing
            break
        
    return