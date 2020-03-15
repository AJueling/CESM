#!/bin/bash
#SBATCH -t 00:25:00
#SBATCH -p normal
#SBATCH -n 19
newgrp imau
cd /home/dijkbio/andre/CESM/src/run
module load pre2019
module load eb
module load Miniconda3
source activate CESM

(
python ../run/run_SST_ctrl_video.py 0 20
)&
(
python ../run/run_SST_ctrl_video.py 20 40
)&
(
python ../run/run_SST_ctrl_video.py 40 60
)&
(
python ../run/run_SST_ctrl_video.py 60 80
)&
(
python ../run/run_SST_ctrl_video.py 80 100
)&
(
python ../run/run_SST_ctrl_video.py 100 120
)&
(
python ../run/run_SST_ctrl_video.py 120 140
)&
(
python ../run/run_SST_ctrl_video.py 140 160
)&
(
python ../run/run_SST_ctrl_video.py 160 180
)&
(
python ../run/run_SST_ctrl_video.py 180 200
)&
(
python ../run/run_SST_ctrl_video.py 200 220
)&
(
python ../run/run_SST_ctrl_video.py 220 240
)&
(
python ../run/run_SST_ctrl_video.py 240 260
)&
(
python ../run/run_SST_ctrl_video.py 260 280
)&
(
python ../run/run_SST_ctrl_video.py 280 300
)&
(
python ../run/run_SST_ctrl_video.py 300 320
)&
(
python ../run/run_SST_ctrl_video.py 320 340
)&
(
python ../run/run_SST_ctrl_video.py 340 360
)&
(
python ../run/run_SST_ctrl_video.py 360 380
)&
wait