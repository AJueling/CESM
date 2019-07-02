#!/bin/bash
#SBATCH -t 20:00:00
#SBATCH -N 1
cd /home/dijkbio/andre/CESM/src
module load eb
module load Miniconda3
source activate CESM
python SST_generation.py had
python SST_generation.py lpd
python SST_generation.py ctrl