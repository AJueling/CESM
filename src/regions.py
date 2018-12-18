# high res ocean regions

import numpy as np
import xarray as xr

from xr_DataArrays import example_file
from paths import file_ex_ocn_rect

ocn_file = example_file('ocn')

# 'ocn' locations
# AMO_area      = {'nlat':slice(1081,1858), 'nlon':slice( 500,1100)}  # North Atlantic (0-60N, 0-80W)
Drake_Passage = {'nlat':slice( 268, 513), 'nlon':410             }
DP_North      = {'nlat':512             , 'nlon':410             }
global_ocean  = {'nlat':slice(   0,2400), 'nlon':slice(   0,3600)}  # global ocean
Nino12        = {'nlat':slice(1081,1181), 'nlon':slice( 200, 300)}  # Niño 1+2 (0-10S, 90W-80W)
Nino34        = {'nlat':slice(1131,1232), 'nlon':slice(3000,3500)}  # Niño 3.4 (5N-5S, 170W-120W)
sinking_area  = {'nlat':slice( 283, 353), 'nlon':slice(1130,1210)}
SOM_area      = {'nlat':slice( 603, 808), 'nlon':slice( 600,1100)}  # (50S-35S, 0E-50W)
# TexT_area     = {'nlat':slice( 427,1858)                         }  # tropic + extratropics (60S-60N)
WGKP_area     = {'nlat':slice(   0, 603), 'nlon':slice( 750,1900)}
WG_center     = {'nlat':321             , 'nlon':877             }  # [64.9S,337.8E]

# 'ocn_rect' locations
# SOM_area_rect = {  'lat': slice(-50,-35),   'lon': slice(0,50)}
gl_ocean_rect = {'t_lat': slice(-80, 90), 't_lon': slice(0,360)}

# 'ocn_low' location
# AMO_area_low  = {'nlat':slice( 187, 353), 'nlon':slice( 284,35)}  # North Atlantic (0-60N, 0-80W)
Nino12_low    = {'nlat':slice( 149, 187), 'nlon':slice( 275, 284)}  # Niño 1+2 (0-10S, 90W-80W)
Nino34_low    = {'nlat':slice( 168, 205), 'nlon':slice( 204, 248)}  # Niño 3.4 (5N-5S, 170W-120W)
# TexT_area_low = {'nlat':slice(  36, 353)                         }  # tropic + extratropics (60S-60N)
SOM_area_low  = {'nlat':slice( 603, 808), 'nlon':slice( 600,1100)}  # 

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
    assert domain in ['ocn', 'ocn_low', 'ocn_rect']
    
    file  = example_file(domain)
    
    if domain in ['ocn', 'ocn_low']:
        RMASK = xr.open_dataset(file, decode_times=False).REGION_MASK
        if mask_nr==0:  # global ocean
            MASK = np.where(RMASK>0, 1, 0)
        else:
            MASK = np.where(RMASK==mask_nr, 1, 0)   
    
    elif domain=='ocn_rect':
        assert mask_nr==0
        RMASK = xr.open_dataset(file, decode_times=False).TEMP[0,:,:]
        MASK = np.where(RMASK>-2, 1, 0)
        
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


def NPacific_mask_rect():
    """ boolean Pacific mask """
    TEMP = xr.open_dataset(file_ex_ocn_rect, decode_times=False).TEMP[0,:,:]
    MASK1 = (TEMP/TEMP).where(TEMP>0, 0)
    MASK1 = MASK1.where(MASK1.t_lat>0  , 0)
    MASK1 = MASK1.where(MASK1.t_lat<66 , 0)
    MASK1 = MASK1.where(MASK1.t_lon>100, 0)
    MASK1 = MASK1.where(MASK1.t_lon<300, 0)
    MASK2 = MASK1.where(MASK1.t_lon>250, 0)
    MASK1 = MASK1.where(MASK1.t_lon<250, 0)
    MASK2 = MASK2.where(2.1*MASK2.t_lat+MASK2.t_lon<300, 0)
    return MASK1+MASK2

def AMO_mask(domain):
    file = example_file(domain)
    TLAT = xr.open_dataset(file, decode_times=False).TLAT
    if domain=='ocn_low':
        MASK = boolean_mask(domain, mask_nr=6)
    elif domain=='ocn':
        MASK = Atlantic_mask(domain)
    MASK = np.where(TLAT> 0, MASK, 0)
    MASK = np.where(TLAT<60, MASK, 0)
    return MASK

def TexT_mask(domain):
    file = example_file(domain)
    TLAT = xr.open_dataset(file, decode_times=False).TLAT
    MASK = boolean_mask(domain, mask_nr=0)
    MASK = np.where(TLAT>-60, MASK, 0)
    MASK = np.where(TLAT< 60, MASK, 0)
    return MASK

def SOM_mask(domain):
    # (50S-35S, 0E-50W)
    file  = example_file(domain)
    TLAT  = xr.open_dataset(file, decode_times=False).TLAT
    TLONG = xr.open_dataset(file, decode_times=False).TLONG
    MASK = boolean_mask(domain, mask_nr=0)
    MASK = np.where(TLAT >-50, MASK, 0)
    MASK = np.where(TLAT <-35, MASK, 0)
    MASK = np.where(TLONG>310, MASK, 0)
    return MASK


def TPI_masks(domain, region_nr):

    file  = example_file(domain)
    TLAT  = xr.open_dataset(file, decode_times=False).TLAT
    TLONG = xr.open_dataset(file, decode_times=False).TLONG
    MASK = boolean_mask(domain, mask_nr=0)
    if region_nr==1:
        # (25N-45N, 140E-145W)
        MASK = np.where(TLAT > 25, MASK, 0)
        MASK = np.where(TLAT < 45, MASK, 0)
        MASK = np.where(TLONG>140, MASK, 0)
        MASK = np.where(TLONG>215, MASK, 0)
    elif region_nr==2:
        # (10S-10N, 170E-90W)
        MASK = np.where(TLAT >-10, MASK, 0)
        MASK = np.where(TLAT < 10, MASK, 0)
        MASK = np.where(TLONG>170, MASK, 0)
        MASK = np.where(TLONG>270, MASK, 0)
    elif region_nr==3:
        # (50S-15S, 150E-160W)
        MASK = np.where(TLAT >-50, MASK, 0)
        MASK = np.where(TLAT <-15, MASK, 0)
        MASK = np.where(TLONG>150, MASK, 0)
        MASK = np.where(TLONG>200, MASK, 0)
    return MASK

def SST_index_bounds(name):
    if name=='TPI1':   boounds = (140,215, 25, 45)
    elif name=='TPI2': boounds = (170,270,-10, 10)
    elif name=='TPI3': boounds = (150,200,-50,-15)
    elif name=='SOM':  boounds = (310,360,-50,-35)
    elif name=='AMO':  boounds = (280,360,  0, 60)
    return bounds