import sys
sys.path.append("..")
from OHC import OHC_integrals
run     = sys.argv[1]
mask_nr = int(sys.argv[2])
domain  = 'ocn'
assert run=='ctrl' or run=='rcp'
assert mask_nr>=0 and mask_nr<13

OHC_integrals(domain=domain, run=run, mask_nr=mask_nr)