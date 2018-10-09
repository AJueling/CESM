import xarray as xr

def xa_vol_int(xa, hpos='TT', vpos='T', pbc=True):
    """volume integral of xarray

    input:
    xa       .. xarray DataArray
    h/vpos   .. horiz./vert. grid position of quantity

    output:
    integral .. float integral
    """
    assert type(xa)==
    if hpos=='TT':
        assert 'TAREA' in xa
    if vpos=='T':
        assert 
    if pbc:
        assert
    
    
    return 0.


def xa_vol_mean(xa):
    """
    
    """
    return


def xa_surf_int(xa, hpos='TT'):
    """surface integral of xarray DataArray
    
    """
    
    return 0.