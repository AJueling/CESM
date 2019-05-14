# code structure

(model output)

(how to deal with grids consistently?)

0. supporting functions and classes
    * IterateCESMoutput

1. generation of files
    (derive_files.py)
    * MakeDerivedFiles
    - yrly averages
    - isolate SST
    - calculate OHC
    * GenerateSSTFields
    - detrended fields

2. analysis of single runs
    (analysis.py)
    * xrAnalaysis
    * FieldAnalysis
    * TimeSeriesAnalysis
    - AMOC
    - spectrum
    - AR(1) process
    
    (analysis_SST.py)
    * IndexAnalysis

3. synthesis of several runs / compare between model runs
    (synthesis.py)
    - common time series
    * TimeSeriesSynthesis ?