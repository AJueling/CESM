# high res ocean regions

import numpy as np
import xarray as xr

from xr_DataArrays import example_file

regions_dict = {-14: 'Caspian_Sea',
                -13: 'Black_Sea',
#                  -1: 'no calculations',
#                   0: 'continents',
                  0: 'Global_Ocean',
                  1: 'Southern_Ocean',
                  2: 'Pacific_Ocean',
                  3: 'Indian_Ocean',
                  4: 'Persian_Gulf',
                  5: 'Red_Sea',
                  6: 'Atlantic_Ocean',
                  7: 'Mediterranean',
                  8: 'Labrador_Sea',
                  9: 'Greenland_Sea',
                 10: 'Arctic_Ocean',
                 11: 'Hudson_Bay',
                 12: 'Baltic_Sea',
               }


def boolean_mask(domain, mask_nr):
    """ 3D boolean xr DataArray """
    assert domain=='ocn'
    
    file  = example_file(domain)
    RMASK = xr.open_dataset(file, decode_times=False).REGION_MASK
    
    if mask_nr==0:  # global ocean
        MASK = np.where(RMASK>0, 1, 0)
    else:
        MASK = np.where(RMASK==mask_nr, 1, 0)
    
    return MASK