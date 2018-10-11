# from src directory call "python -m unittest test.test_xr_integrate"
# https://stackoverflow.com/questions/1896918/running-unittest-with-typical-test-directory-structure

import unittest
import numpy as np
import xarray as xr
import xr_integrate
from xr_integrate import xr_vol_int , xr_vol_mean , xr_int_along_axis
from xr_integrate import xr_surf_int, xr_surf_mean

from xr_DataArrays import create_xr_DataArray, generate_xr_DZ, generate_xr_AREA
from paths import file_ex_ocn_ctrl, file_ex_ocn_rect


class TestXrIntegrate(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        # runs once before all test
        # set up mock xarray data array
        # 1. simplest version
        z = np.array([0,1])
        y = np.array([0,1,2])
        x = np.array([0,1,2,3])
        cls.A = xr.DataArray(data=np.ones((2,3,4)),
                              coords={'z': z, 'y': y, 'x': x},
                              dims=('z', 'y', 'x') )
        cls.B = xr.DataArray(data=np.ones((3,4)),
                              coords={'y': y, 'x': x},
                              dims=('y', 'x') )

        # 2. full sized version, simple entries
            # C     .. 3D xr Dataset of high res CESM output
            # D     .. 3D xr DataArray of ones
            # DZT   .. 3D xr DataArray of depths in [m]
            # TAREA .. 2D xr DataArray of T-cell surface ares [m^2]
            
        # 2.1 hires ocean fields
        dim_names = ('nlat','nlon','z_t')
        cls.D, cls.C, imt, jmt, km = create_xr_DataArray(file_ex_ocn_ctrl,
                                                         dim_names,
                                                         n=3,
                                                         fill=1)
        cls.DZT = generate_xr_DZ(case='ocn_hires_fbc')
        assert cls.DZT.sum().item()==18962679618.534042
        # [m]; sum of all depths: not physically sensible
        cls.TAREA = generate_xr_AREA(case='ocn_hires')
        assert cls.TAREA.sum().item()==360713803999209.06  # [m^2]

        # 2.2 atm fields

        
    
    def test_xr_vol_int(self):
        # simple arrays
        A, B = self.A, self.B
        lats = np.arange(-89.5,90,1)

        self.assertEqual(xr_vol_int(A, B, A, levels=False, zonal=False),
                         24.)
        
        self.assertEqual(xr_vol_int(A, B, A, levels=False, zonal=True )[0],
                         24.)
        self.assertTrue(np.all(xr_vol_int(A, B, A, levels=False, zonal=True )[1].coords['y_bins']==lats))
        
        self.assertEqual(xr_vol_int(A, B, A, levels=True , zonal=False)[0],
                         24.)
        self.assertTrue(np.all(xr_vol_int(A, B, A, levels=True , zonal=False)[1]==12.))
        
        self.assertEqual(xr_vol_int(A, B, A, levels=True , zonal=True )[0],
                         24.)
        self.assertTrue(np.all(xr_vol_int(A, B, A, levels=True , zonal=True )[1]==12.))
        self.assertTrue(np.any(xr_vol_int(A, B, A, levels=True, zonal=True)[2]==4.))
        self.assertTrue(np.any(xr_vol_int(A, B, A, levels=True, zonal=True)[3]==1.))
        
        
        # ocean arrays
        self.assertAlmostEqual(xr_vol_int(self.D, self.TAREA, self.DZT),\
                               1.3660833e+18,
                               delta=1e11)  # vol of the ocean

    def test_xr_int_along_axis(self):
        A, B = self.A, self.B
        self.assertTrue(np.all(xr_int_along_axis(A, A, 0)==2.))
        self.assertTrue(np.all(xr_int_along_axis(A, A, 1)==3.))
        self.assertTrue(np.all(xr_int_along_axis(A, A, 2)==4.))
        
        self.assertTrue(np.all(xr_int_along_axis(B, B, 0)==3.))
        self.assertTrue(np.all(xr_int_along_axis(B, B, 1)==4.))
        
#     def test_xr_vol_int_regional(self):
#         self.assertEqual(xr_integrate.xr_vol_int(self.A, self.B, self.A),
#                          xr_integrate.xr_vol_int_regional(self.A, self.B, self.A, self.)) 


#     def test_find_regional_coord_extent(self):
#         # create boolean array
#         pass


    def test_xa_vol_mean(self):
        A, B = self.A, self.B
        self.assertEqual(xr_vol_mean(A, B, A), 1.)
        
        
    def test_xr_surf_int(self):
        # simple arrays
        B = self.B
        self.assertEqual(xr_surf_int(B, B), 12.)
        
        # ocean arrays
        self.assertAlmostEqual(xr_surf_int(self.D[0,:,:], self.TAREA),\
                               3.6071380e+14,
                               delta=1e7)  # surface of the ocean
 

    def test_xr_surf_mean(self):
        B = self.B
        self.assertEqual(xr_surf_mean(B, B), 1.)


#     def test_zonal_int(self):


#     def test_xr_zonal_mean(self):
        # zonal means os 1 field should be 1
#         pass

        
if __name__ == '__main__':
    # if this is not included, the script needs to be executed via
    # 'python unittest -m test_calc'
    unittest.main()
    
    
# extra notes:

# there two ways to check whether the ValueError gets raised
#self.assertRaises(ValueError, calc.divide, 10, 0)
#with self.assertRaises(ValueError):  # this is called a context manager
#    calc.divide(10,0)