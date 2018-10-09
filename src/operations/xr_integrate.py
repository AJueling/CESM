import numpy as np
import xarray as xr

from analysis.paths import file_ex_ocn_hires, file_ex_ocn_rect



def create_xr_DataArray(file, n=3, fill=0):
    """creates xr DataArray object with n last coords/dims of file
    
    input:
    file .. contains xr Dataset with specific dims/coords to be copied
    n    .. number of dimensions: 2 -> horizontal 2D array, 3 -> full 3D array
    
    output:
    DA   .. n-dimensional xr DataArray object with zeros
    """
    C = xr.open_dataset(file, decode_times=False)
    imt,jmt,km = C.dims['nlon'], C.dims['nlat'], C.dims['z_t']
    
    if n==2:
        if fill==0:  array = np.zeros((jmt,imt))
        if fill==1:  array = np.ones((jmt,imt))
        DA = xr.DataArray(data=array,
                         coords={'nlat': C.coords['nlat'], 'nlon': C.coords['nlon']},
                         dims=('nlat', 'nlon') )
        
    if n==3:
        if fill==0:  array = np.zeros((km,jmt,imt))
        if fill==1:  array = np.ones((km,jmt,imt))
        DA = xr.DataArray(data=array,
                         coords={'z_t': C.coords['z_t'], 'nlat': C.coords['nlat'], 'nlon': C.coords['nlon']},
                         dims=('z_t', 'nlat', 'nlon') )
    
    return DA, C, imt, jmt, km



def generate_xr_DZ(hpos='TT', vpos='T', pbc=False):
    """ builds 3D xr DataArray of cell depths in [m]
    
    input:
    hpos .. horizontal position
    vpos .. vertical position
    pbc  .. partial bottom cells
    
    output:
    DZT  .. 3D xr DataArray object with depths in [m]
    """
    DZ, C, imt, jmt, km = create_xr_DataArray(file_ex_ocn_hires, 3)

    for k in range(km):
        if pbc==False and vpos=='T' and hpos=='TT':
            DZ[k,:,:] = np.where(C.KMT[:,:]<=k, DZ[k,:,:], C.dz[k]/100)

    return DZ



def generate_xr_AREA(hpos='TT', k=0):
    """ builds 2D xr DataArray of surface area [m^2]
    
    input:
    hpos .. horizontal position on grid
    k    .. level at which surface is to be created
    
    output:
    AREA .. 2D xr DataArray
    """
    AREA, C, imt, jmt, km = create_xr_DataArray(file_ex_ocn_hires, 2)
    if hpos=='TT':
        AREA[:,:] = np.where(C.KMT[:,:]<=k, AREA[:,:], C.TAREA/1e4)
    
    return AREA



def xr_vol_int(xa, AREA, DZ):
    """volume integral of xarray

    input:
    xa       .. xr DataArray
    AREA     .. 2D xr DataArray
    DZ       .. 3D xr DataArray 

    output:
    integral .. float integral
    """
    assert type(xa)==xr.core.dataarray.DataArray
    assert type(AREA)==xr.core.dataarray.DataArray
    assert type(DZ)==xr.core.dataarray.DataArray
    assert np.shape(AREA)==np.shape(xa)[-2:]
    assert np.shape(DZ)==np.shape(xa)[-3:]
    
    integral = np.sum(xa*AREA*DZ)
    
    return integral.item()



def xr_vol_mean(xa, AREA, DZ):
    """
    (to be written)
    """
    return



def xr_surf_int(xa, AREA):
    """surface integral of xarray DataArray
    
    input:
    xa       .. 2D xr DataArray
    AREA     .. 2D xr DataArray of surface area
    
    output:
    integral .. float integrated
    """
    assert type(xa)==xr.core.dataarray.DataArray
    assert type(AREA)==xr.core.dataarray.DataArray
    assert np.shape(xa)==np.shape(AREA)
    assert len(np.shape(xa))==2
    
    integral = np.sum(xa*AREA)
    
    return integral.item()