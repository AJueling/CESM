#!/bin/bash
#SBATCH -t 1:00:00
#SBATCH -N 1
cd /home/dijkbio/andre/CESM/src/regrid
module load eb
module load Miniconda3
source activate CESM
python regrid_yrly_TEMP_PD.py