""" copies versions of specific notebooks stripped of their output to the analysis folder """

import os
import sys
import shutil
import io
from nbstripout import strip_output, read, write, NO_CONVERT

destination = '../'
assert os.path.exists(destination)

# copying stripped .ipynb files
files = ['currents',
         'data_overview',
         'geometry',
         'GMST',
         'OHC',
         'SHF',
         'SST',
         'timeseries',
         'winds',
        ]

for x in files:
    fn = f'{x}.ipynb'
    with io.open(fn, 'r', encoding='utf8') as f:
        nb = read(f, as_version=NO_CONVERT)
    nb = strip_output(nb, keep_output=False, keep_count=False)
    with io.open(destination+fn, 'w', encoding='utf8') as f:
        write(nb, f)
    print(f'copied {fn}')
