import os
import numpy as np
import xarray as xr

from paths import file_ex_ocn_hires, file_ex_ocn_rect, file_ex_atm_hires
from constants import R_earth

def create_xr_DataArray(file, dim_names, n=3, fill=0):
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
    assert os.path.exists(file)==True
    assert len(dim_names)==3
    assert n==2 or n==3
    assert fill==0 or fill==1
    
    (lat, lon, z) = dim_names
    C = xr.open_dataset(file, decode_times=False)
    imt, jmt, km = C.dims[lon], C.dims[lat], C.dims[z]
    
    if n==2:
        if fill==0:  array = np.zeros((jmt,imt))
        if fill==1:  array = np.ones((jmt,imt))
        coords = {lat: C.coords[lat], lon: C.coords[lon]}
        DA = xr.DataArray(data=array,
                          coords=coords,
                          dims=(lat, lon) )
        
    if n==3:
        if fill==0:  array = np.zeros((km,jmt,imt))
        if fill==1:  array = np.ones((km,jmt,imt))
        coords = {z: C.coords[z], lat: C.coords[lat], lon: C.coords[lon]}
        DA = xr.DataArray(data=array,
                          coords=coords,
                          dims=(z, lat, lon) )
    
    return DA, C, imt, jmt, km



def generate_xr_DZ(case):
    """ builds 3D xr DataArray of cell depths in [m]
    
    input:
    hpos .. horizontal position
    vpos .. vertical position
    pbc  .. partial bottom cells
    
    output:
    DZT  .. 3D xr DataArray object with depths in [m]
    """
    if case=='ocn_hires_fbc' or case=='ocn_hires_pbc':
        dim_names = ('nlat','nlon','z_t')
        file = file_ex_ocn_hires
        fill = 0
    elif case=='atm':
        dim_names = ('lat','lon','lev')
        file = file_ex_atm_hires
        fill = 1
    else:
        raise ValueError('case argument not known')
    
    DZ, C, imt, jmt, km = create_xr_DataArray(file=file,
                                              dim_names=dim_names,
                                              n=3,
                                              fill=0)
    
    if case=='ocn_hires_fbc':
        for k in range(km):
            DZ[k,:,:] = np.where(C.KMT[:,:]<=k, DZ[k,:,:], C.dz[k]/100)
    elif case=='ocn_hires_pbc':
        print('not implemented yet')
    elif case=='atm':
        print('how is this implemented?')
        
    return DZ



def generate_xr_AREA(case):
    """ builds 2D xr DataArray of surface area [m^2]
    
    input:
    hpos .. horizontal position on grid
    k    .. level at which surface is to be created
    
    output:
    AREA .. 2D xr DataArray [m^2]
    """
    if case=='ocn_hires':
        dim_names = ('nlat','nlon','z_t')
        file = file_ex_ocn_hires
        k=0  # area at surface
    elif case=='atm':
        dim_names = ('lat','lon','lev')
        file = file_ex_atm_hires
    else:
        raise ValueError('case argument not known')
    AREA, C, imt, jmt, km = create_xr_DataArray(file=file,
                                                dim_names=dim_names,
                                                n=2,
                                                fill=0)
    
    if case=='ocn_hires':
        AREA[:,:] = np.where(C.KMT[:,:]<=k, AREA[:,:], C.TAREA/1e4)
    elif case=='atm':
        dy = C.lat[1].item()-C.lat[0].item()
        print(dy)
        nx, ny = len(C.lon), len(C.lat)
        lat_N = (C.lat[0]+dy/2)*np.pi/180
        AREA[0 ,:] = spher_surf_element(R_earth, 2*np.pi/nx, lat_N, -np.pi/2)
        lat_S = (C.lat[-1]-dy/2)*np.pi/180
        AREA[-1,:] = spher_surf_element(R_earth, 2*np.pi/nx, np.pi/2, lat_S)     
        for j in range(1, ny-1):
            lat_S = (C.lat[j]-dy/2)*np.pi/180
            lat_N = (C.lat[j]+dy/2)*np.pi/180
            AREA[j,:] = spher_surf_element(R_earth, 2*np.pi/nx, lat_N, lat_S)
        assert np.isclose(np.max(AREA.values),
                          R_earth**2 * 2*np.pi/nx * np.pi/ny, rtol=1e-2)
    
    return AREA

    
def spher_surf_element(r, dtheta, lat_N, lat_S):
    return r**2 * dtheta * (np.sin(lat_N)-np.sin(lat_S))