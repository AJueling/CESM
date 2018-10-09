# from src directory call "python -m unittest test.test_xarray_integrate"
# https://stackoverflow.com/questions/1896918/running-unittest-with-typical-test-directory-structure

import unittest
import xarray as xr
from analysis import paths
from operations import xarray_integrate

from paths import file_ex_ocn_hires, file_ex_ocn_hires_rect


class TestXarrayIntegrate(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        # runs once before all test
        
        # set up mock xarray data array
        # 1. simplest version
        
        # 2. full sized version, simple entries
        
        # 3. 
        
        print('running setUpClass')
        
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
#         # pass
    
#     def tearDown(self):
#         # runs after every single test
#         # could be used to delete files that are created
#         pass
    
    def test_xa_vol_int(self):
        self.assertEqual(xarray_integrate.xa_vol_int(3), 0.)
#         self.assertEqual(calc.add(-1,1), 0)
#         self.assertEqual(calc.add(-1,-1), -2)
#         print('running test_add')
        # testing edge cases
        # then copying them for all functions
        
    def test_xa_vol_mean(self):
        
    def test_xa_surf_int(self):
#         self.assertEqual(calc.subtract(10,5), 5)
#         self.assertEqual(calc.subtract(-1,1), -2)
#         self.assertEqual(calc.subtract(-1,-1), 0)

#     def test_multiply(self):
#         self.assertEqual(calc.multiply(10,5), 50)
#         self.assertEqual(calc.multiply(-1,1), -1)
#         self.assertEqual(calc.multiply(-1,-1), 1)
        
#     def test_divide(self):
#         self.assertEqual(calc.divide(10,5), 2)
#         self.assertEqual(calc.divide(-1,1), -1)
#         self.assertEqual(calc.divide(-1,-1), 1)
#         self.assertEqual(calc.divide(5,2), 2.5)
        # added last test in case simple division was replaced by floor division
        
        # two ways to check whether the ValueError gets raised
        #self.assertRaises(ValueError, calc.divide, 10, 0)
#         with self.assertRaises(ValueError):
#             # this is called a context manager
#             calc.divide(10,0)
        
if __name__ == '__main__':
    # if this is not included, the script needs to be executed via
    # 'python unittest -m test_calc'
    unittest.main()