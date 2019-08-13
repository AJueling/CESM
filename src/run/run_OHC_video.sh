#!/bin/bash
#SBATCH -t 02:10:00
#SBATCH -p normal
#SBATCH -N 1
cd /home/dijkbio/andre/CESM/src/run
module load eb
module load Miniconda3
source activate CESM
python run_OHC_video.py