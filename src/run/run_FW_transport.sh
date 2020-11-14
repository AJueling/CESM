#!/bin/bash
#SBATCH -t 03:30:00
#SBATCH -p fat  # 256GB RAM
#SBATCH -N 2 --ntasks-per-node=2
newgrp imau
cd /home/dijkbio/andre/CESM/src/run
module load pre2019
module load eb
module load Miniconda3
source activate CESM

(
python ../FW_transport.py lpd 154 601
python ../FW_transport.py ctrl 1 101
)&
(
python ../FW_transport.py ctrl 101 201
)&
(
python ../FW_transport.py ctrl 201 301
)&
(
python ../FW_transport.py lr1 2000 2101
python ../FW_transport.py rcp 2000 2101
)&

wait