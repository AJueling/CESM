#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
cd /home/dijkbio/andre/CESM/src
module load eb
module load Miniconda3
source activate CESM
python SST_regression.py lpd
python SST_regression.py ctrl
