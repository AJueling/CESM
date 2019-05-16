# code structure

(model output)

(how to deal with grids consistently?)



a. generation of files
    (derive_files.py)
    * MakeDerivedFiles
    - yrly averages
    - isolate SST
    - calculate OHC
    * GenerateSSTFields
    - detrended fields

b. analysis of single runs
    (analysis.py)
    * xrAnalysis
    * FieldAnalysis
    * TimeSeriesAnalysis
    - AMOC
    - spectrum
    - AR(1) process
    
    (analysis_SST.py)
    * IndexAnalysis

c. synthesis of several runs / compare between model runs
    (synthesis.py)
    - common time series
    * TimeSeriesSynthesis ?
    

d. shared supporting functions and classes
    * IterateCESMoutput
    - CESM file name
    - paths
    
e. specific functions for particular analyses