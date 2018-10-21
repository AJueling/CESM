#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH -N 9 --ntasks-per-node=3
#SBATCH -p normal

cd /home/dijkbio/andre/CESM/src/run

module load eb
module load Miniconda3
source activate CESM

(python run_OHC_integrals.py ctrl  0) &
(python run_OHC_integrals.py ctrl  1) &
(python run_OHC_integrals.py ctrl  2) &
(python run_OHC_integrals.py ctrl  3) &
(python run_OHC_integrals.py ctrl  4) &
(python run_OHC_integrals.py ctrl  5) &
(python run_OHC_integrals.py ctrl  6) &
(python run_OHC_integrals.py ctrl  7) &
(python run_OHC_integrals.py ctrl  8) &
(python run_OHC_integrals.py ctrl  9) &
(python run_OHC_integrals.py ctrl 10) &
(python run_OHC_integrals.py ctrl 11) &
(python run_OHC_integrals.py ctrl 12) &

(python run_OHC_integrals.py rcp  0) &
(python run_OHC_integrals.py rcp  1) &
(python run_OHC_integrals.py rcp  2) &
(python run_OHC_integrals.py rcp  3) &
(python run_OHC_integrals.py rcp  4) &
(python run_OHC_integrals.py rcp  5) &
(python run_OHC_integrals.py rcp  6) &
(python run_OHC_integrals.py rcp  7) &
(python run_OHC_integrals.py rcp  8) &
(python run_OHC_integrals.py rcp  9) &
(python run_OHC_integrals.py rcp 10) &
(python run_OHC_integrals.py rcp 11) &
(python run_OHC_integrals.py rcp 12) &

wait