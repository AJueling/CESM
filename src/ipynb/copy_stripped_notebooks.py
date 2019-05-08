""" copies versions of specific notebooks stripped of their output to the analysis folder """

import os
import sys
import shutil
import io
from nbstripout import strip_output, read, write, NO_CONVERT

destination = '../notebooks/'
assert os.path.exists(destination)

# copying stripped .ipynb files
files = ['AHC',
         'BSF',
         'CICE',
         'CURRENTS',
         'data_overview',
         'FORCING',
         'geometry',
         'GMST_obs_ensembles',
         'GMST',
         'HEAT',
         'HIATUS',
         'MOC',
         'MXL',
         'OHC-lowres',
         'OHC',
         'OSNAP',
         'paper',
         'Poster',
         'REGR',
         'regrid',
         'SHF',
         'SST',
         'SST_AMO',
         'SST_ENSO',
         'SST_GMST_regression',
         'SST_obs',
         'SST_PDO',
         'SST_regression_maps',
         'SST_SOM',
         'SST_TPI',
         'testing_significance',
         'timeseries',
         'TOA',
         'WIND',
        ]

for x in files:
    fn = f'{x}.ipynb'
    with io.open(fn, 'r', encoding='utf8') as f:
        nb = read(f, as_version=NO_CONVERT)
    nb = strip_output(nb, keep_output=False, keep_count=False)
    with io.open(destination+fn, 'w', encoding='utf8') as f:
        write(nb, f)
    print(f'copied {fn}')
