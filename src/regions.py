# high res ocean regions

import numpy as np
import xarray as xr

from xr_DataArrays import example_file

SOM_area_rect = {'lat': slice(-50,-35), 'lon': slice(0,50)}

WGKP_area = {'nlat': slice(0,603), 'nlon': slice(750,1900)}

sinking_area = {'nlat': slice(250,350), 'nlon': slice(1100,1200)}

Drake_Passage = {'nlat':slice(268,513), 'nlon':410}

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


def combine_mask(MASK, numbers):
    """
    combines submasks of MASK into single boolean mask
    
    input:
    MASK    .. xr DataArray with integer numbers for masks
    numbers .. list of IDs of submasks to be combined
    """
    assert len(numbers)>1
    
    NEW_MASK = xr.where(MASK==numbers[0], 1, 0)
    for number in numbers[1:]:
        NEW_MASK = xr.where(MASK==number, 1, NEW_MASK)
    
    return NEW_MASK


def Atlantic_mask(domain):
    assert domain=='ocn'
    
    ds = xr.open_dataset(example_file(domain), decode_times=False)
    ATLANTIC_MASK = combine_mask(ds.REGION_MASK, [6,8,9])
    
    return ATLANTIC_MASK