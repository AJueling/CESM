# Multidecadal Variability paper - data pipeline

## Fig. 1:  global annual SHF
1. annual average SHF files with `timeseries.ipynb`
2. combine files into single file with `AnalyzeBudget().surface_heat_flux(run=run)`

## Fig. 2:  SST index time series

1. create annual SST files
2. `DS().generate_monthly_SST_files(run='lc1')`
2. apply `SST_index_creation.py` checks existence and generates intermediate files if not present
    1. generates 
    
## Fig. 3:  index regression maps
1. also generated in `SST_index_creation.py`

## Fig. 4:  SST index spectra
1. also generated in `SST_index_creation.py`

## Fig. 5:  GMST/TOA-SHF spectra
1. annual average atmospheric files
2. execute function in `GMST.ipynb`

## Fig. 6:  SHF spectra
1. created with SHF integrals of Fig. 1

## Fig. 7:  OHC(t,lat)
1. TEMP_PD annual average files in `timeseries.ipynb`
2. `ac_derivation_OHC.py`

## Fig. 8:  OHC(t,z)
(created with Fig. 7 data)
