import os
import sys
import dask
import numpy as np
import scipy as sp
import xarray as xr

sys.path.append("..")
from tqdm import tqdm_notebook
from LFCA import CESM_xlfca
from paths import path_prace, path_data, path_results
from filters import lowpass, Lanczos, deseasonalize
from regions import boolean_mask, mask_box_in_region
from eofs.xarray import Eof
from xr_DataArrays import xr_AREA
from ab_derivation_SST import DeriveSST as DS
from bb_analysis_timeseries import AnalyzeTimeSeries as ATS

run   = sys.argv[1]
basin = sys.argv[2]
dsdt  = sys.argv[3]

CESM_xlfca(run, basin, dsdt)