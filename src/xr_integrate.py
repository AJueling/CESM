import numpy as np
import xarray as xr


def xr_vol_int(xa, AREA, DZ):
    """ volume integral of xarray

    input:
    xa         .. 3D xr DataArray with data to be integrated
    AREA       .. 2D xr DataArray of cell areas
    DZ         .. 3D xr DataArray of cell depths

    output:
    int_levels .. integrals of each level
    integral   .. float integral
    """
    assert type(xa)==xr.core.dataarray.DataArray
    assert len(np.shape(xa))==3
    assert type(AREA)==xr.core.dataarray.DataArray
    assert type(DZ)==xr.core.dataarray.DataArray
    assert np.shape(AREA)==np.shape(xa)[-2:]
    assert np.shape(DZ)==np.shape(xa)[-3:]
    
    km = len(xa[:,0,0])
    int_levels = np.zeros((km))
    for k in range(km):
        int_levels[k] = np.sum(xa[k,:,:]*AREA[:,:]*DZ[k,:,:]).item()
    integral = np.sum(int_levels)
    
    return integral, int_levels



def xr_vol_int_regional(xa, AREA, DZ, MASK):
    """ volumen integral with  regional MASK
    
    input:
    xa, AREA, DZ         .. same as in 'xr_vol_int'
    MASK                 .. 2D xr DataArray of booleans with the same dimensions as xa 
    
    output:
    integral, int_levels .. same as in 'xr_vol_int'
    
    """
    
    assert type(xa)==xr.core.dataarray.DataArray
    assert type(AREA)==xr.core.dataarray.DataArray
    assert type(DZ)==xr.core.dataarray.DataArray
    assert np.shape(AREA)==np.shape(xa)[-2:]
    assert np.shape(DZ)==np.shape(xa)[-3:]
    assert np.dtype(MASK)==np.dtype('bool')

    # determine min/max i/j of masked region
    (imin, imax, jmin, jmax) = find_regional_coord_extent(MASK)
    
    xa_reg   = xa.where(MASK)[:,jmin:jmax+1,imin:imax+1]
    AREA_reg = AREA.where(MASK)[jmin:jmax+1,imin:imax+1]
    DZ_reg   = DZ.where(MASK)[:,jmin:jmax+1,imin:imax+1]
    
    integral, int_levels = xr_vol_int(xa_reg, AREA_reg, DZ_reg)
   
    return integral, int_levels



def find_regional_coord_extent(MASK):
    """ finds coordinates of a boolean mask
    
    input:
    MASK .. 2D xr DataArray of booleans
    
    output:
    (imin, imax, jmin, jmax) .. lon/lat extent of True area
    """
    assert type(MASK)==xr.core.dataarray.DataArray
    
    jmin = np.where(MASK)[0].min()
    jmax = np.where(MASK)[0].max()
    imin = np.where(MASK)[1].min()
    imax = np.where(MASK)[1].max()
    
    return (imin, imax, jmin, jmax)



def xr_vol_mean(xa, AREA, DZ):
    """ 
    (to be written)
    """
    return



def xr_surf_int(xa, AREA):
    """ surface integral of xarray DataArray
    
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



def xr_zonal_int(xa, AREA, dx, lat_name):
    """ integral over dx wide latitude bins
        
    input:
    xa          .. 2D xr DataArray
    AREA        .. 2D xr DataArray
    dx          .. width of latitude bands
    lat_name    .. xa/AREA coordinate name of the latitude variable
    
    output:
    xa_zonal_int  .. 1D xr DataArray
    lat_center  .. np.ndarray with bin centers
    """
    
    assert type(xa)==xr.core.dataarray.DataArray
    assert type(AREA)==xr.core.dataarray.DataArray
    assert len(np.shape(xa))==2
    assert np.shape(xa)==np.shape(AREA)
    assert dx>(xa[lat_name][-1]-xa[lat_name][0])/len(xa[:,0])
    
    lat_bins = np.arange(-90, 90+dx, dx)
    lat_center = np.arange(-90+dx/2, 90, dx)
    
    xa_new = xa*AREA
    xa_zonal_int = xa_new.groupby_bins(lat_name, lat_bins, labels=lat_center).sum()
    
    return xa_zonal_int, lat_center



def xr_zonal_mean(xa, AREA, dx, lat_name):
    """ area weighted mean over dx wide latitude bins
        
    input:
    xa          .. 2D xr DataArray
    AREA        .. 2D xr DataArray
    dx          .. width of latitude bands
    lat_name    .. xa/AREA coordinate name of the latitude variable
    
    output:
    xa_zonal_mean .. 1D xr DataArray
    lat_center  .. np.ndarray with bin centers
    """
    
    assert type(xa)==xr.core.dataarray.DataArray
    assert type(AREA)==xr.core.dataarray.DataArray
    assert len(np.shape(xa))==2
    assert np.shape(xa)==np.shape(AREA)
    assert dx>180/len(AREA[0,:])
    
    xa_zonal_int, lat_center = xr_zonal_int(xa, AREA, dx, lat_name)
    AREA_zonal_int, lat_center = xr_zonal_int(AREA/AREA, AREA, dx, lat_name)
    
    xa_zonal_mean = xa_zonal_int/AREA_zonal_int

    return xa_zonal_mean, lat_center