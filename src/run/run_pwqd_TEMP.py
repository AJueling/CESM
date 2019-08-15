import os
import sys
import datetime
sys.path.append("..")

# import dask
# from dask.distributed import Client, LocalCluster
from aa_derivation_fields import DeriveField as DF

run = sys.argv[1]
# run = 'lpd'
import dask
# from dask.distributed import Client, LocalCluster
# cluster = LocalCluster(n_workers=8)
# client = Client(cluster)

if __name__=="__main__":
    assert run in ['ctrl', 'lpd']
#     cluster = LocalCluster(n_workers=1)
#     client = Client(cluster)
    
    print(f'running pwqd TEMP run={run}')
    print(f'{datetime.datetime.now()}\n\n')
    DF(run).make_pwqd_TEMP_files()
    print(f'\n\nfinished at\n{datetime.datetime.now()}')