# TargeTox - Biological Network Target-based Drug Toxicity Risk Prediction

!!! 27/10/2018 Update currently in progress and the new version is due to be released shortly !!!

## Required libraries:

  * “catboost” – machine learning library used for implementation of the main model
  * “sets” - used for efficient Gene Ontology indexing in implementation of Functional Impact score
  
## Recommended libraries 

  * “STRINGdb” - as the method only accepts STRING database identifiers, this library can be used to convert other types for use with TargeTox

## Usage instructions

Download and extract the script and supporting data in the "R" sub-directory. Inside R, load the script by running:

`source('TargeTox.r')`

The main method can be run with a call to TargeTox function which takes the following arguments:
 (1) a vector of String ids of protein targets 
 (2-4) route of administration values (administration.oral, administration.parenteral, administration.topical) where
  values are 0 - not administered via that route; 1 - administered and 2 - unclassified; all defined in equivalent way to ChEMBL database.
 (5-6) Lower and upper plasma protein binding values in a 0-100 range (can be missing)

 The method returns a single TargeTox score, with higher score indicating higher toxicity risk. Note that these returned values are not probabilities and will not necessarily be within 0-1 range.
 
 The file also contains code for two supporting functions that return shortest diffusion state distance to reference nodes (minDsd) and Functional Impact score (functionalImpact). Both of these functions take a vector of STRING ids as input.
