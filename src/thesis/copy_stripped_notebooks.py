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
files = ['2_forcings',
         '2_simulations',
         '6.1_indicators',
         '6.2_ECS',
         '6.3_OHC',
         '6.3_SAT',
         '6.4_radiation',
         '6.5_seaice',
        ]

for x in files:
    fn = f'{x}.ipynb'
    with io.open(fn, 'r', encoding='utf8') as f:
        nb = nbstripout.read(f, as_version=nbstripout.NO_CONVERT)
    nb = nbstripout.strip_output(nb, keep_output=False, keep_count=False)
    with io.open(destination+fn, 'w', encoding='utf8') as f:
        nbstripout.write(nb, f)
    print(f'copied {fn}')
