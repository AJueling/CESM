import os
import numpy as np
import xarray as xr

from paths import file_ex_ocn_hires, file_ex_ocn_rect


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
        DA = xr.DataArray(data=array,
                          coords={lat: C.coords[lat], lon: C.coords[lon]},
                          dims=(lat, lon) )
        
    if n==3:
        if fill==0:  array = np.zeros((km,jmt,imt))
        if fill==1:  array = np.ones((km,jmt,imt))
        DA = xr.DataArray(data=array,
                          coords={z: C.coords[z], lat: C.coords[lat], lon: C.coords[lon]},
                          dims=(z, lat, lon) )
    
    return DA, C, imt, jmt, km



def generate_xr_DZ(file, hpos='TT', vpos='T', pbc=False):
    """ builds 3D xr DataArray of cell depths in [m]
    
    input:
    hpos .. horizontal position
    vpos .. vertical position
    pbc  .. partial bottom cells
    
    output:
    DZT  .. 3D xr DataArray object with depths in [m]
    """
    dim_names = ('nlat','nlon','z_t')
    DZ, C, imt, jmt, km = create_xr_DataArray(file=file,
                                              dim_names=dim_names,
                                              n=3,
                                              fill=0)
    
    assert 'dz' in C
    
    for k in range(km):
        if pbc==False and vpos=='T' and hpos=='TT':
            DZ[k,:,:] = np.where(C.KMT[:,:]<=k, DZ[k,:,:], C.dz[k]/100)

    return DZ



def generate_xr_AREA(file, hpos='TT', k=0):
    """ builds 2D xr DataArray of surface area [m^2]
    
    input:
    hpos .. horizontal position on grid
    k    .. level at which surface is to be created
    
    output:
    AREA .. 2D xr DataArray [m^2]
    """
    dim_names = ('nlat','nlon','z_t')
    AREA, C, imt, jmt, km = create_xr_DataArray(file=file,
                                                dim_names=dim_names,
                                                n=2,
                                                fill=0)
    if hpos=='TT':
        AREA[:,:] = np.where(C.KMT[:,:]<=k, AREA[:,:], C.TAREA/1e4)
    
    return AREA