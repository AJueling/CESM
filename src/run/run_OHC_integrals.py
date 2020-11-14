import sys
import datetime
sys.path.append("..")

from ac_derivation_OHC import DeriveOHC as DO

run  = sys.argv[1]
year = int(sys.argv[2])

assert run in ['ctrl', 'rcp', 'hq', 'lpd', 'lpi', 'lc1', 'lq']

print(f'running OHC_integrals run={run} year={year}')
print(f'{datetime.datetime.now()}\n\n')
DO().generate_OHC_files(run=run, year=year)
print(f'\n\nfinished at\n{datetime.datetime.now()}')