#!/bin/bash
#SBATCH -t 3:30:00
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH -p normal

cd /home/dijkbio/andre/CESM/src/run

module load eb
module load Miniconda3
source activate CESM

(python run_yrly_avg.py 1) &
(python run_yrly_avg.py 2) &
(python run_yrly_avg.py 3) &
(python run_yrly_avg.py 4) &

wait