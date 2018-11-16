import sys
import datetime
sys.path.append("..")

from OHC import OHC_integrals

run     = sys.argv[1]
mask_nr = int(sys.argv[2])
domain  = 'ocn'

assert run=='ctrl' or run=='rcp'
assert mask_nr>=0 and mask_nr<13

print(f'running OHC_integrals run={run} region={mask_nr}')
print(f'{datetime.datetime.now()}\n\n')
OHC_integrals(domain=domain, run=run, mask_nr=mask_nr)
print(f'\n\nfinished at\n{datetime.datetime.now()}')