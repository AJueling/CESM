import sys
import datetime
sys.path.append("..")

from aa_derivation_fields import DeriveField as DF

run = sys.argv[1]

if __name__=="__main__":
    assert run in ['ctrl', 'lpd']
    print(f'running pwqd TEMP run={run}')
    print(f'{datetime.datetime.now()}\n\n')
    DF(run).make_pwqd_TEMP_files()
    print(f'\n\nfinished at\n{datetime.datetime.now()}')