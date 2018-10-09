import os
import sys
import shutil
import io
from nbstripout import strip_output, read, write, NO_CONVERT

destination = '../analysis/'
assert os.path.exists(destination)

# copying stripped .ipynb files
files = ['data_overview',
        ]
for x in files:
    fn = f'{x}.ipynb'
    with io.open(fn, 'r', encoding='utf8') as f:
        nb = read(f, as_version=NO_CONVERT)
    nb = strip_output(nb, keep_output=False, keep_count=False)
    with io.open(destination+fn, 'w', encoding='utf8') as f:
        write(nb, f)
    print(f'copied {fn}')