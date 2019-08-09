import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import mtspec
from pandas.tools.plotting import autocorrelation_plot
from itertools import combinations

from maps import regr_map
from paths import path_samoc, path_results
from constants import A_earth
from analysis import TimeSeriesAnalysis, FieldAnalysis
# from xr_regression import lag_linregress_3D

    
class FieldSynthesis(FieldAnalysis):
    """ compares xr fields """
    def __init__(self, run):
        self.run = run
        
#     def spatial_
