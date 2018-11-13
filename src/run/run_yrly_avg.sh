#!/bin/bash
#SBATCH -t 4:00:00
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH -p normal

(
cd /home/dijkbio/andre/CESM/src/run
module load eb
module load Miniconda3
source activate CESM
python run_yrly_avg.py 4
)&
wait