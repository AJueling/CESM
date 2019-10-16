#!/bin/bash
#SBATCH -t 06:00:00
#SBATCH -p normal
#SBATCH -N 1
cd /home/dijkbio/andre/CESM/src/run
newgrp imau
module load pre2019
module load eb
module load Miniconda3
source activate CESM

(
python run_inter_basin_fluxes.py ctrl OHC
)&
(
python run_inter_basin_fluxes.py ctrl SALT
)&
wait