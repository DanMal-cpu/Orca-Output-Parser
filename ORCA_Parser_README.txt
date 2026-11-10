
ORCA Data Parser - Main Script
===============================

This script provides the main interface for parsing and analyzing ORCA output files.
It handles both standard calculations and TD-DFT calculations with Multiwfn integration.

Example ORCA output filenames: Fe3+_Mo_Keggin_H-t-butyl-Peroxo.out

All you need to do is provide your file of choice, and the parser will detect the type of
calculation used and act as such:

Single Point, Optimization, or Frequency 
- properties and parameters are collected and analysed
- individual parameters can be printed by calling via the console (see full list below) 
- save in a pair of txt files, one containing all data, the other a more focused file with key
findings.

TD-DFT
- generates a UV-Vis spectra, with an adjustable sigma value
- finds high intensity and local peaks, retrieves NTO data for these peaks
- generates relevant NTO hole, NTO electron, and CDD cube files for all key transitions
- saves all generated data in neat, ordered tranche of folders
