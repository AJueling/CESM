#!/bin/bash
#SBATCH -t 05:00:00
#SBATCH -p normal
#SBATCH -N 1
cd /home/dijkbio/andre/CESM/src/run
module load pre2019
module load eb
module load Miniconda3
source activate CESM
python run_OHC_integrals.py ctrl  # 45 mins
python run_OHC_integrals.py lpd   # 26 mins