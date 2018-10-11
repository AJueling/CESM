import numpy as np
import xarray as xr


def xr_vol_int(xa, AREA, DZ, levels=False, zonal=False):
    """ volume integral of xarray *[m^3]

    input:
    xa                  .. 3D xr DataArray with data to be integrated
    AREA                .. 2D xr DataArray of cell areas
    DZ                  .. 3D xr DataArray of cell depths
    levels              .. option to output results for all level
    zonal               .. option to output zonal integrals

    output:
    integral            .. float integral
    int_levels          .. integrals of each level
    xa_zonal_int        .. 1D array of vert.+zonally integrated quantity
    xa_zonal_level_int  .. 2D (km, lat_bin) *[m^2] (integrated in depth and lon)
    xa_zonal_level_mean .. 2D (km, lat_bin) *[m^1] 
                           (weighted by bottom cell depth)
    """
    assert type(xa)==xr.core.dataarray.DataArray
    assert len(np.shape(xa))==3
    assert type(AREA)==xr.core.dataarray.DataArray
    assert type(DZ)==xr.core.dataarray.DataArray
    assert np.shape(AREA)==np.shape(xa)[-2:]
    assert np.shape(DZ)==np.shape(xa)[-3:]
    
    if zonal==True:
        dx = 1  # latitude bin width
        if np.shape(DZ)==(2,3,4):           # simple test case
            lat_name = 'y'
        elif np.shape(DZ)==(42,2400,3600):  # hires ocean
            lat_name = 'nlat'
        elif np.shape(DZ)==(30,384,576):    # atm fields
            lat_name = 'lat'
        else:
            raise ValueError('unknown shape: lat_name not implemented')
        assert lat_name in DZ.coords
    
    if levels==False:
        integral = np.sum(xa[:,:,:]*AREA[:,:]*DZ[:,:,:]).item()
        
        if zonal==False:  # just global integral
            return integral
        
        elif zonal==True:
            xa_vert = xr_int_along_axis(xa, DZ, 0)
            xa_zonal_int = xr_zonal_int(xa_vert, AREA, dx, lat_name)
            return integral, xa_zonal_int
        
    elif levels==True:
        km = len(xa[:,0,0])
        int_levels = np.zeros((km))
        for k in range(km):
            int_levels[k] = np.sum(xa[k,:,:]*AREA[:,:]*DZ[k,:,:]).item()
        integral = np.sum(int_levels)
        
        if zonal==False:
            return integral, int_levels
        
        if zonal==True:
            ONES = AREA.copy()
            ONES[:,:] = 1.
            for k in range(km):
                xa_zonal_int = xr_zonal_int(xa[k,:,:]*DZ[k,:,:], AREA, dx, lat_name)
                DZ_zonal_int = xr_zonal_int(DZ[k,:,:]          , ONES, dx, lat_name)
                if k==0:
                    xa_zonal_level_int = np.zeros((km, len(xa_zonal_int)))
                    xa_zonal_level_mean = np.zeros((km, len(xa_zonal_int)))
                xa_zonal_level_int[k,:] = xa_zonal_int
                xa_zonal_level_mean[k,:] = xa_zonal_int/DZ_zonal_int
            return integral, int_levels, xa_zonal_level_int, xa_zonal_level_mean

        
        
def xr_int_along_axis(xa, DZ, axis):
    """ integral of xr DataArray along a specific axis 
    
    input:
    xa   .. 3D xr DataArray of quantity to be integrated
    DZ   .. 3D xr DataArray of vertical cell extents [m]
    axis .. int axis to be integrated over
    
    output:
    int  .. 2D xr DataArray of integrated quantitity
    """
    assert type(axis)==np.dtype(int)
    assert np.shape(xa)==np.shape(DZ)
    assert axis<=len(np.shape(xa))
    
    integral = np.sum(xa*DZ, axis=axis)
    
    return integral
        
    

def xr_vol_int_regional(xa, AREA, DZ, MASK):
    """ volumen integral with regional MASK
    
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
    """ mean over quantity stored in xa
    
    input:
    xa   .. 3D xr DataArray of quantity
    AREA .. 2D xr DataArray of cell surface area
    DZ   .. 3D xr DataArray of cell depths
    
    output:
    mean .. (float)
    """
    assert type(xa)==xr.core.dataarray.DataArray
    assert type(DZ)==xr.core.dataarray.DataArray
    assert type(AREA)==xr.core.dataarray.DataArray
    assert np.shape(xa)==np.shape(DZ)
    assert np.shape(xa[0,:,:])==np.shape(AREA)    
    
    integral    = xr_vol_int(xa, AREA, DZ, levels=False, zonal=False)
    ONES        = xa.copy()
    ONES[:,:,:] = 1.
    volume      = xr_vol_int(ONES, AREA, DZ, levels=False, zonal=False)
    mean        = integral/volume
    
    return mean



def xr_surf_int(xa, AREA):
    """ surface integral of xarray DataArray *[m^2]
    
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



def xr_surf_mean(xa, AREA):
    """ mean over a surface *[1] 
    
    input:
    xa   .. 2D xr DataArray of quantity
    AREA .. 2D xr DataArrayof cell surfaces
    
    output:
    mean .. (float) mean of quantity in xa
    """
    assert type(xa)==xr.core.dataarray.DataArray
    assert type(AREA)==xr.core.dataarray.DataArray
    assert np.shape(xa)==np.shape(AREA)
    assert len(np.shape(xa))==2

    integral  = xr_surf_int(xa, AREA)
    ONES      = xa.copy()
    ONES[:,:] = 1.
    surface   = xr_surf_int(ONES, AREA)
    mean      = integral/surface
    
    return mean



def xr_zonal_int(xa, AREA, dx, lat_name):
    """ integral over dx wide latitude bins
        
    input:
    xa          .. 2D xr DataArray
    AREA        .. 2D xr DataArray
    dx          .. width of latitude bands
    lat_name    .. xa/AREA coordinate name of the latitude variable
    
    output:
    xa_zonal_int  .. 1D xr DataArray
    
    lat centers can be accessed through xa_zonal_int.coords[f'{lat_name}_bins']
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
    
    return xa_zonal_int



def xr_zonal_mean(xa, AREA, dx, lat_name):
    """ area weighted mean over dx wide latitude bins
        
    input:
    xa            .. 2D xr DataArray
    AREA          .. 2D xr DataArray
    dx            .. width of latitude bands
    lat_name      .. xa/AREA coordinate name of the latitude variable
    
    output:
    xa_zonal_mean .. 1D xr DataArray
    """
    
    assert type(xa)==xr.core.dataarray.DataArray
    assert type(AREA)==xr.core.dataarray.DataArray
    assert len(np.shape(xa))==2
    assert np.shape(xa)==np.shape(AREA)
    assert dx>180/len(AREA[0,:])
    
    xa_zonal_int = xr_zonal_int(xa, AREA, dx, lat_name)
    AREA_zonal_int = xr_zonal_int(AREA/AREA, AREA, dx, lat_name)
    
    xa_zonal_mean = xa_zonal_int/AREA_zonal_int

    return xa_zonal_mean