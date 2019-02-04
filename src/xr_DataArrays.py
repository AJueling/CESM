import os
import numpy as np
import xarray as xr

from paths import file_ex_ocn_ctrl, file_ex_ocn_rect, file_ex_atm_ctrl,\
                  file_ex_ice_rcp, file_ex_atm_lpd, file_ex_ocn_lpd, file_ex_atm_lpi,\
                  path_samoc, path_results, file_geometry, file_HadISST
from constants import R_earth
from read_binary import read_binary_2D_double



def create_xr_DataArray(domain, n=3, fill=0):
    """creates xr DataArray object with n last coords/dims of file
    
    input:
    file       .. contains xr Dataset with specific dims/coords to be copied
    names      .. tuple of lat, lon and height/depth dimension names
    n          .. number of dimensions: 2 -> horizontal 2D array, 3 -> full 3D array
    fill       .. returned array with 0's or 1's
    
    output:
    DA         .. n-dimensional xr DataArray object with zeros or ones
    C          .. xr Dataset
    imt,jmt,km .. dimensions of lon,lat,depth/height
    """
    
    assert n in [2, 3]
    assert fill in [0, 1]
    
    (z, lat, lon) = depth_lat_lon_names(domain)
    file = example_file(domain)
    
    C = xr.open_dataset(file, decode_times=False)
    imt, jmt = C.dims[lon], C.dims[lat]
    if domain in ['ocn', 'ocn_low', 'ocn_rect', 'atm', 'atm_f19', 'atm_f09', 'ice']:
        km = C.dims[z]
    elif domain=='ocn_had':
        km = None
    
    if n==2:
        if fill==0:  array = np.zeros((jmt,imt))
        if fill==1:  array = np.ones((jmt,imt))
        coords = {lat: C.coords[lat], lon: C.coords[lon]}
        da = xr.DataArray(data=array,
                          coords=coords,
                          dims=(lat, lon) )
        
    if n==3:
        if fill==0:  array = np.zeros((km,jmt,imt))
        if fill==1:  array = np.ones((km,jmt,imt))
        coords = {z: C.coords[z], lat: C.coords[lat], lon: C.coords[lon]}
        da = xr.DataArray(data=array,
                          coords=coords,
                          dims=(z, lat, lon) )
    
    return da, C, imt, jmt, km



def xr_DZ(domain, grid='T'):
    """ builds 3D xr DataArray of cell depths in [m]
    
    input:
    domain .. (str) 'ocn_hires_fbc'
    
    output:
    DZT    .. 3D xr DataArray object with depths in [m]
    """
    assert domain in ['ocn', 'ocn_low', 'ocn_rect']
    
    if grid=='U':  file_path = f'{path_samoc}/geometry/DZU_{domain}.nc'
    else:          file_path = f'{path_samoc}/geometry/DZT_{domain}.nc'
        
    if os.path.exists(file_path):
        DZ = xr.open_dataarray(file_path)
        
    else:
        DZ, C, imt, jmt, km = create_xr_DataArray(domain=domain, n=3, fill=0)
            
        if domain=='ocn':  # partial bottom cells
            PBC = read_binary_2D_double(file_geometry, 3600, 2400, 1)  # [lon, lat]
            
            if grid=='U':
                for k in range(km):
                    DZ[k,:,:] = np.where(C.KMU[:,:]>k , C.dz[k]/100   , DZ[k,:,:])
                    # technically there should be the 
            else:  # T-grid
                for k in range(km):
                    DZ[k,:,:] = np.where(C.KMT[:,:]>k , C.dz[k]/100   , DZ[k,:,:])
                    DZ[k,:,:] = np.where(C.KMT[:,:]==k, PBC[:,:].T/100, DZ[k,:,:])
                                
        elif domain=='ocn_low':
            for k in range(km):
                DZ[k,:,:] = np.where(C.PD[0,k,:,:]>0, (C.z_w_bot[k]-C.z_w_top[k])/1e2, 0)
            
        elif domain=='ocn_rect':
            for k in range(km):
                DZ[k,:,:] = np.where(C.PD[k,:,:]>0, C.w_dep[k+1]-C.w_dep[k], DZ[k,:,:])
        
        DZ.to_netcdf(file_path)
        
    return DZ



def xr_AREA(domain):
    """ builds 2D xr DataArray of surface area [m^2]
    
    input:
    hpos .. horizontal position on grid
    k    .. level at which surface is to be created
    
    output:
    AREA .. 2D xr DataArray [m^2]
    """
    assert domain in ['ocn', 'ocn_low', 'ocn_rect', 'ocn_had', 'atm', 'atm_f19', 'atm_f09', 'ice']
    
    AREA, C, imt, jmt, km = create_xr_DataArray(domain=domain, n=2, fill=0)
    
    if domain in ['ocn', 'ocn_low', 'ice']:  # TAREA of cells are written out
        AREA[:,:] = C.TAREA/1e4
        
    elif domain in ['ocn_had', 'atm', 'atm_f19', 'atm_f09']:  # regular, rectangular grids
        (z, lat, lon) = depth_lat_lon_names(domain)
        dy = abs(C[lat][1].item()-C[lat][0].item())
        nx, ny = len(C[lon]), len(C[lat])
        jmin, jmax = 1, ny-1
        lat_N = (-90+dy/2)*np.pi/180
        AREA[0 ,:] = spher_surf_element(R_earth, 2*np.pi/nx, lat_N, -np.pi/2)
        lat_S = (90-dy/2)*np.pi/180
        AREA[-1,:] = spher_surf_element(R_earth, 2*np.pi/nx, np.pi/2, lat_S)
            
        for j in range(jmin, jmax):
            lat_S = (C[lat][j]-dy/2)*np.pi/180
            lat_N = (C[lat][j]+dy/2)*np.pi/180
            AREA[j,:] = spher_surf_element(R_earth, 2*np.pi/nx, lat_N, lat_S)

        # ensure calculated area max is within 1 order of magnitude of a naive cell area
#         assert np.isclose(np.log10(np.max(AREA.values)),\
#                           np.log10(R_earth**2 * 2*np.pi/nx * np.pi/ny), rtol=1)
                
        if domain=='ocn_had':  # lat elements saved North to South
            AREA = abs(AREA)
            
    elif domain=='ocn_rect':  # rect grid with different lat.-diffs.
        ds = xr.open_dataset(file_ex_ocn_rect, decode_times=False)
        nx = len(ds.t_lon)
        for j, l in enumerate(ds.t_lat[:-1]):
            AREA[j,:] = spher_surf_element(R_earth, 2*np.pi/nx, ds.t_lat[j+1]*np.pi/180, ds.t_lat[j]*np.pi/180)

    return AREA



def xr_HTN(domain):
    """ zonal length of grid cells
    
    returns
    HTN .. 2D xr DataArray
    """
    assert domain in ['ocn', 'ocn_low', 'ocn_rect', 'atm'] 
    
    (z, lat, lon) = depth_lat_lon_names(domain)
    HTN, C, imt, jmt, km = create_xr_DataArray(domain=domain, n=2, fill=0)
    
    if domain in ['ocn', 'ocn_low']:
        # this is not exactly the zonal length of the T-cell at its center
        # however, the error introduced is smaller than 1%
        # also it is not directed zonally in the Northern NH due to hte tripolar grid
        HTN = C.HTN
        
    elif domain=='ocn_rect':
        n_lat, n_lon = len(HTN.lat), len(HTN.lon)
        for j in range(n_lat):
            HTN[j,:] = zonal_length(HTN.lat[j].item(), n_lon)

    return HTN


def xr_LATS(domain):
    """ latitudes of grid cells
    
    returns
    LATS .. 2D xr DataArray
    """
    assert domain in ['ocn', 'ocn_low']#, 'ocn_rect', 'atm'] 
    
    (z, lat, lon) = depth_lat_lon_names(domain)
    
    if domain in ['ocn', 'ocn_low']:
        LATS = xr.open_dataset(file_ex_ocn_ctrl, decode_times=False).TLAT

#     elif domain=='ocn_rect':
#         HTN, C, imt, jmt, km = create_xr_DataArray(domain=domain, n=2, fill=0)
#         n_lat, n_lon = len(HTN.lat), len(HTN.lon)
#         for j in range(n_lat):
#             HTN[j,:] = zonal_length(HTN.lat[j].item(), n_lon)

    return LATS



def depth_lat_lon_names(domain):
    """ old syntax """
    return dll_dims_names(domain)


def dll_dims_names(domain):
    """ names for different model domain dimensions 
    
    input:
    domain .. (str)
    
    output:
    dll    .. tuple of strings of depth, lat, lon dimension names
    """
    assert domain in ['ocn', 'ocn_low', 'ocn_rect', 'ocn_had', 'atm', 'ice', 'atm_f19', 'atm_f09']
    
    if domain in ['ocn', 'ocn_low']:
        dll = ('z_t', 'nlat', 'nlon')
    elif domain=='ocn_rect':
        dll = ('depth_t', 't_lat', 't_lon')
    elif domain=='ocn_had':
        dll = (None, 'latitude', 'longitude')
    elif domain in ['atm', 'atm_f19', 'atm_f09']:
        dll = ('lev', 'lat', 'lon')    
    return dll


def dll_coords_names(domain):
    """ names for different model domain coordinates 
    
    input:
    domain .. (str)
    
    output:
    dll    .. tuple of strings of depth, lat, lon coordinate names
    """
    assert domain in ['ocn', 'ocn_low', 'ocn_rect', 'ocn_had', 'atm', 'ice', 'atm_f19', 'atm_f09']
    
    if domain=='ocn':
        dll = ('z_t', 'TLAT', 'TLONG')
    if domain=='ocn_low':
        dll = ('z_t', 'nlat', 'nlon')
    elif domain=='ocn_rect':
        dll = ('depth_t', 't_lat', 't_lon')
    elif domain=='ocn_had':
        dll = (None, 'latitude', 'longitude')
    elif domain in ['atm', 'atm_f19', 'atm_f09']:
        dll = ('lev', 'lat', 'lon')
    
    return dll



def dll_from_arb_da(da):
    """ finds dimension names from arbitrary xr DataArray via its size"""
    assert type(da)==xr.core.dataarray.DataArray
    shape = np.shape(da)[-1]
    assert shape in [3600, 900, 320]
    
    if shape==3600:
        dll = depth_lat_lon_names('ocn')
    elif shape==900:
        dll = depth_lat_lon_names('ocn_rect')
    elif shape==320:
        dll = depth_lat_lon_names('ocn_low')
    elif shape==360:
        dll = depth_lat_lon_names('ocn_had')
    return dll



def example_file(domain):
    """ example of output file for a given domain """
    assert domain in ['ocn', 'ocn_low', 'ocn_rect', 'ocn_had', 'atm', 'ice', 'atm_f19', 'atm_f09']

    if   domain=='ocn':       file = file_ex_ocn_ctrl
    elif domain=='ocn_low':   file = file_ex_ocn_lpd
    elif domain=='ocn_rect':  file = file_ex_ocn_rect
    elif domain=='ocn_had':   file = file_HadISST
    elif domain=='atm':       file = file_ex_atm_ctrl
    elif domain=='atm_f19':   file = file_ex_atm_lpi
    elif domain=='atm_f09':   file = file_ex_atm_lpd
    elif domain=='ice':       file = file_ex_ice_rcp
    
    assert os.path.exists(file)==True
    
    return file

    
    
def spher_surf_element(r, dtheta, lat_N, lat_S):
    """ surface area of element of sphere """
    return r**2 * dtheta * (np.sin(lat_N)-np.sin(lat_S))



def zonal_length(lat, nlon):
    """ length of zonal 1/nlon segment at latitude lat"""
    return R_earth * 2*np.pi/nlon * np.cos(lat*np.pi/180)



def xr_DXU(domain):
    """ U-cell zonal length [m] """
    assert domain=='ocn'
    fn = f'{path_samoc}/geometry/DXU.nc'
    
    if os.path.exists(fn)==True:
        da = xr.open_dataarray(fn)
    else:
        f = example_file(domain)
        da = xr.open_dataset(f, decode_times=False).DXU/1e2
        da.to_netcdf(fn)
    
    return da


def xr_DYU(domain):
    """ U-cell zonal length [m] """
    assert domain=='ocn'
    fn = f'{path_samoc}/geometry/DYU.nc'
    
    if os.path.exists(fn)==True:
        da = xr.open_dataarray(fn)
    else:
        f = example_file(domain)
        da = xr.open_dataset(f, decode_times=False).DYU/1e2
        da.to_netcdf(fn)
    
    return da