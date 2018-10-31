#!/bin/bash
#SBATCH -t 8:00:00
#SBATCH -N 1 --ntasks-per-node=2
#SBATCH -p normal

(
cd /home/dijkbio/andre/CESM/src/run
module load eb
module load Miniconda3
source activate CESM
python run_yrly_avg.py  9
)
&
(
cd /home/dijkbio/andre/CESM/src/run
module load eb
module load Miniconda3
source activate CESM
python run_yrly_avg.py 10
)
&
wait