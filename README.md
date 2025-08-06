#  Biphasic generalised modelling and Bayesian analysis framework for coral reef recovery analysis

Tools for Modeling and inference for reef recovery modelling across the Great Barrier Reef. Modelling is motivated by results of the paper,

DJ Warne, KA Crossman, GEM Heron, JA Sharp, W Jin, PP-Y Wu, MJ Simpson, K Mengersen, J-C Ortiz (2024) Mathematical modelling and uncertainty quantification for analysis of biphasic coral reef recovery patterns. *ArXiv.org*. https://arxiv.org/abs/2406.19591

## Summary

This repository contains tools to perform the following:

1. Extraction, filtering and processing of data from the LTMP and MMP for use in the modelling framework;
2. Tools for direct identification of two-phase recovery based on two-segment regression;
3. Definition of various ODE-based models for coral reef recovery. This account for various features such as two-phase recovery, multispecies interactions, and uncertainty in the fundamental logistic growth form of coral growth; 
4. Bayesian inference methods for calibration and uncertainty quantification of generalised recovery models;
5. Prective sampling (both prior predictive and posterior predictive) for model testing
6. Example processing pipelines.


## Developers

 The following developers implemented these tools:
 
 - David J. Warne [1,2,3] (david.warne@qut.edu.au)
 - Grace E. M. Heron [1,3] (g.heron@qut.edu.au)
 
 Affilations:
 
   1. School of Mathematical Sciences, Faculty of Science, Queensland University of Technology 
   2. Centre for Data Science, Queensland University of Technology
   3. ARC Centre of Excellence for Mathematical and Statistical Frontiers

## Acknowledgements

This software is part of a joint project that has been funded by:

1. ARC Centre of Excellence for Mathematical and Statistical Frontiers (ACEMS)
2. Centre for Data Science, Queensland University of Technology
3. Australian Institute for Marine Science (AIMS)



## Requirements

**Software Requirements:** 

 - R Version 4.2.3 (2023-3-13) Stortstop Beagle
 - RStudio 2023.03.0+386 "Cherry Blossom" Release (3c53477afb13ab959aeb5b34df1f10c237b256c3, 2023-03-09) for Windows

**Important R packages:**

   - tidyverse
   - deSolve 
   - mcmc 
   - adaptMCMC 
   - coda 
   - doFuture 


**Note**

Larger analysis runs will likely require HPC resources. We utilised QUT's HPC cluster which is administered by the eResearch Office at QUT.

## Contents

This repository contains the following:

- `LTMPTools/` Main library of functions for data processing, modelling, inference and prediction
- `models/` Model definition files implementing various biphasic Richards' recovery models
- `data_processing/` Example scripts for data extraction and filtering 
- `pipelines/` Contains example pipeline demonstrating usage of the framework for inference, predition and plotting

