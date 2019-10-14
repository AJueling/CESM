import sys
import datetime
sys.path.append("..")

from be_analysis_budget import AnalyzeBudget

run        = sys.argv[1]
quantitity = sys.argv[2]

assert run in ['ctrl', 'lpd']
assert quantity in ['SALT', 'OHC']

print(f'running OHC_integrals run={run}')
print(f'{datetime.datetime.now()}\n\n')
AnalyzeBudget().all_transports(run=run, quantity=q)
print(f'\n\nfinished at\n{datetime.datetime.now()}')