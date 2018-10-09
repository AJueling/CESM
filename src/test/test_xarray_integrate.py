# from src directory call "python -m unittest test.test_xarray_integrate"
# https://stackoverflow.com/questions/1896918/running-unittest-with-typical-test-directory-structure

import unittest
import numpy as np
import xarray as xr
from operations import xr_integrate
from analysis.paths import file_ex_ocn_hires, file_ex_ocn_rect


class TestXarrayIntegrate(unittest.TestCase):
    
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
        cls.D, cls.C, imt, jmt, km = xr_integrate.create_xr_DataArray(file_ex_ocn_hires, n=3, fill=1)
        
        cls.DZT = xr_integrate.generate_xr_DZ(hpos='TT', vpos='T', pbc=False)
        assert cls.DZT.sum().item()==18962679618.534042  # [m]; sum of all depths: not physically sensible
        
        cls.TAREA = xr_integrate.generate_xr_AREA(hpos='TT', k=0)
        assert cls.TAREA.sum().item()==360713803999209.06  # [m^2]

#     @classmethod
#     def tearDownClass(cls):
#         pass
    
#     def setUp(self):
#         # this setup method will be run before every single test
#         # can create toy instances of a class, for example,
#         # which can then be used by several tests
#         # variables/class instanced need to be set as TestClass instances
#         # eg. self.emp1 = Employee('John', 'Doe')
#         # where Employee creates a class instance
#         print('running setUp')


#     def tearDown(self):
#         # runs after every single test
#         # could be used to delete files that are created
#         pass
    
    def test_xa_vol_int(self):
        self.assertEqual(xr_integrate.xr_vol_int(self.A, self.B, self.A), 24.)
        self.assertEqual(xr_integrate.xr_vol_int(self.D, self.TAREA, self.DZT),\
                         1.366083344045523e+18)  # vol of the ocean
        # testing edge cases
        # then copying them for all functions
        
#     def test_xa_vol_mean(self):
        
    def test_xa_surf_int(self):
        self.assertEqual(xr_integrate.xr_surf_int(self.B, self.B), 12.)
        self.assertEqual(xr_integrate.xr_surf_int(self.D[0,:,:], self.TAREA),\
                         360713803999209.06)  # surface of the ocean
        
        # two ways to check whether the ValueError gets raised
        #self.assertRaises(ValueError, calc.divide, 10, 0)
#         with self.assertRaises(ValueError):
#             # this is called a context manager
#             calc.divide(10,0)
        
if __name__ == '__main__':
    # if this is not included, the script needs to be executed via
    # 'python unittest -m test_calc'
    unittest.main()