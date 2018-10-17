# DEPENDENCIES

All computations are done on Cartesius of Surfsara.

The python 3.7 virtual environment used to run the analysis is saved as [] in the documentation folder.

The 'run' scripts that enable scheduled (parallel) computing on Cartesius use sbatch.


# WORKFLOW

1. time averaging
    selected fields of monthly output .nc files to yearly averaged .nc files
2. analysis scripts



# FILE OVERVIEW

## /doc: Documentation
(non-tracked files, including model documentation pdf's)


## /src/  (source code)

python files

- constants.py
- grid.py
- maps.py
- OHC.py
- paths.py
- plotting.py
- read_binary.py
- regions.py
- timesries.py
- xr_DataArrays
- xr_integration
- xr_regression

ipython notebooks
- currents
- data_overview
- geometry
- GMST
- OHC
- SHF
- SST
- timeseries
- winds


### /src/ipynb/  (ipython notebooks including output)

- copy_stripped_notebooks
    makes a copy of the ipynb's with output and strips it out

### /src/run/  (scripts and python files to run in parallel)

- 
- 

### /src/test/  (test suite)

- test_xr_integrals.py

