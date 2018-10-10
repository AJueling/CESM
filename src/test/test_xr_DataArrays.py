# from src directory call "python -m unittest test.test_xr_DataArrays"
# https://stackoverflow.com/questions/1896918/running-unittest-with-typical-test-directory-structure

import unittest
import numpy as np
import xarray as xr
import xr_DataArrays

from paths import file_ex_ocn_hires, file_ex_ocn_rect


class TestXrDataArrays(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        pass
    
    def test_create_xr_DataArray(self):
        pass
    
    def test_generate_xr_DZ(self):
        pass
    
    def test_generate_xr_AREA(self):
        pass

if __name__ == '__main__':
    # if this is not included, the script needs to be executed via
    # 'python unittest -m test_calc'
    unittest.main()