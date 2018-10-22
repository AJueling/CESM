#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH -N 1 --ntasks-per-node=1
#SBATCH -p normal

cd /home/dijkbio/andre/CESM/src/run

module load eb
module load Miniconda3
source activate CESM


(python run_OHC_integrals.py rcp 11) &

wait
