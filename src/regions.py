# high res ocean regions

import numpy as np
import xarray as xr

from xr_DataArrays import example_file

ocn_file = example_file('ocn')

# 'ocn' locations
Drake_Passage = {'nlat':slice( 268, 513), 'nlon':410             }
DP_North      = {'nlat':512             , 'nlon':410             }
global_ocean  = {'nlat':slice(   0,2400), 'nlon':slice(   0,3600)}  # global ocean
Nino12        = {'nlat':slice(1081,1181), 'nlon':slice( 200, 300)}  # Niño 1+2 (0-10S, 90W-80W)
Nino34        = {'nlat':slice(1131,1232), 'nlon':slice(3000,3500)}  # Niño 3.4 (5N-5S, 170W-120W)
sinking_area  = {'nlat':slice( 283, 353), 'nlon':slice(1130,1210)}
SOM_area      = {'nlat':slice( 603, 808), 'nlon':slice( 600,1100)}  # (50S-35S, 0E-50W)
WGKP_area     = {'nlat':slice(   0, 603), 'nlon':slice( 750,1900)}
WG_center     = {'nlat':321             , 'nlon':877             }  # [64.9S,337.8E]

# 'ocn_rect' locations
SOM_area_rect = {'lat': slice(-50,-35), 'lon': slice(0,50)}

# 'atm' locations
Uwind_eq_Pa   = {'lat':slice(-6,6), 'lon':slice(180,200)}


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
