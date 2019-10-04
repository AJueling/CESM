#!/bin/bash

# activates the conda CESM virtual environment and creates a jupyter kernel that can be accessed remotely
module load pre2019
module load eb
module load Miniconda3
source activate CESM
echo 'to access jupyter server remotely: ssh -N -f -L localhost:8890:localhost:8891 dijkbio@cartesius.surfsara.nl'
echo 'in browser type: localhost:8890'
jupyter lab --no-browser --port=8891

