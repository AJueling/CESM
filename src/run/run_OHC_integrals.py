import sys
import datetime
sys.path.append("..")

from ac_derivation_OHC import DeriveOHC as DO

run     = sys.argv[1]

assert run in ['ctrl', 'rcp', 'lpd', 'lpi']

print(f'running OHC_integrals run={run}')
print(f'{datetime.datetime.now()}\n\n')
DO().generate_OHC_files(run=run)
print(f'\n\nfinished at\n{datetime.datetime.now()}')