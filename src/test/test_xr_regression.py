# from src directory call "python -m unittest test.test_xr_regression"
# https://stackoverflow.com/questions/1896918/running-unittest-with-typical-test-directory-structure

import unittest
import numpy as np
import xarray as xr
from xr_regression import xr_quadtrend


# self.assertEqual('foo'.upper(), 'FOO')
# self.assertTrue('FOO'.isupper())
# self.assertFalse('Foo'.isupper())


class TestXrRegression(unittest.TestCase):
    
    # @classmethod
    def setUp(self):
        a = np.arange(5)+1
        b = (a**2).astype(float)
        self.A = xr.DataArray(data=b,
                        dims=('time'),
                        coords={'time':a},
                        )
        self.B = xr.DataArray(data=np.column_stack([b]*5),
                        dims=('time', 'x'),
                        coords={'time':a, 'x':a},
                        )
        self.C = xr.DataArray(data=np.column_stack([b]*10).reshape(5,5,2),
                        dims=('time', 'x', 'y'),
                        coords={'time':a, 'x':a, 'y':[1,2]},
                        )
        pass
    
    def test_xr_quadtrend(self):
        # throw error if no time coordinate not first
        def are_close(x, y):
            self.assertIsNone(xr.testing.assert_allclose(x,y))

        are_close(self.A, self.A)
        are_close(xr_quadtrend(self.A), self.A)
        are_close(xr_quadtrend(self.B), self.B)

        print(self.C)
        print(xr_quadtrend(self.C))

        are_close(xr_quadtrend(self.C), self.C)
        pass

if __name__ == '__main__':
    # if this is not included, the script needs to be executed via
    # 'python unittest -m test_calc'
    unittest.main()