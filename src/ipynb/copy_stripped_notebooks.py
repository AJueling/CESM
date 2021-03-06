""" copies versions of specific notebooks stripped of their output to the analysis folder """

import os
import sys
import shutil
import io
import nbstripout
# from nbstripout import strip_output, read, write, NO_CONVERT

destination = '../notebooks/'
assert os.path.exists(destination)

# copying stripped .ipynb files
files = ['AHC',
         'BSF',
         'CICE',
         'CURRENTS',
         'data_overview',
         'EN4',
         'filters',
         'FLUX_OHC_SALT',
         'FORCING',
         'FW_budget',
         'FW_fluxes',
         'FW_fluxes_BS_Med',
         'geometry',
         'GMST_obs_ensembles',
         'GMST',
         'HEAT',
         'HIATUS',
         'MOC',
         'MXL',
         'OHC-ctrl_vs_lpd',
         'OHC-detrending',
         'OHC-lowres',
         'OHC-observations',
         'OHC-phasing',
         'OHC-videos',
         'OHC',
         'OSNAP',
         'paper',
         'paper3',
         'paper_EPJ',
         'Poster',
         'presentation_BBOS19',
         'REGR',
         'regrid_0.1_to_0.4',
         'regrid_tutorial',
         'SFWF',
         'SHF',
         'SPECTRA',
         'SST',
         'SST_AMO',
         'SST_detrending',
         'SST_ENSO',
         'SST_GMST_regression',
         'SST_GMST_spectra',
         'SST_indices',
         'SST_obs',
         'SST_PDO',
         'SST_regression',
         'SST_SOM',
         'SST_TPI',
         'testing_MSSA',
         'testing_SSA',
         'testing_significance',
         'timeseries',
         'TOA',
         'WIND',
         'Woosok',
        ]

for x in files:
    fn = f'{x}.ipynb'
    with io.open(fn, 'r', encoding='utf8') as f:
        nb = nbstripout.read(f, as_version=nbstripout.NO_CONVERT)
    nb = nbstripout.strip_output(nb, keep_output=False, keep_count=False)
    with io.open(destination+fn, 'w', encoding='utf8') as f:
        nbstripout.write(nb, f)
    print(f'copied {fn}')
