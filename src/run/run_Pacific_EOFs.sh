#!/bin/bash
#SBATCH -t 01:00:00
#SBATCH -p normal
#SBATCH -n 9
newgrp imau
cd /home/dijkbio/andre/CESM/src/run
module load eb
module load Miniconda3
source activate CESM

(
python run_Pacific_EOFs.py ctrl 38S
)&
(
python run_Pacific_EOFs.py ctrl Eq
)&
(
python run_Pacific_EOFs.py ctrl 20N
)&
(
python run_Pacific_EOFs.py lpd 38S
)&
(
python run_Pacific_EOFs.py lpd Eq
)&
(
python run_Pacific_EOFs.py lpd 20N
)&
(
python run_Pacific_EOFs.py had 38S
)&
(
python run_Pacific_EOFs.py had Eq
)&
(
python run_Pacific_EOFs.py had 20N
)&
wait



# for run in ['ctrl', 'lpd', 'had']:
#     for extent in ['38S', 'Eq', '20N']:
#         print('(')
#         print(f'python run_Pacific_EOFs.py {run} {extent}')
#         print(')&')