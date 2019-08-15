#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -p fat
#SBATCH -N 1
cd /home/dijkbio/andre/CESM/src/run
module load eb
module load Miniconda3
source activate CESM
python run_pwqd_TEMP.py ctrl
# python run_pwqd_TEMP.py lpd