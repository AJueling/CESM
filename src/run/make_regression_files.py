import sys
sys.path.append("..")
import scipy as sp
import numpy as np
import xarray as xr

from synthesis import IndexAnalysis

for index in ['AMO']:#['TPI', 'SOM', 'AMO']:
    print(index)
    IA = IndexAnalysis(index)
    IA.make_regression_files()
