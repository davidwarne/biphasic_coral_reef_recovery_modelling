##
# @file script_main_init_analysis_pipeline.R
# @brief complete initial pipeline for reproducibility.
#
# @author David J. Warne (david.warne@qut.edu.au)
#           ARC Centre of Excellence for Mathematical and Statistical Frontiers
#           School of Mathematical Sciences 
#           Faculty of Science
#           Queensland University of Technology
# 
#

# setup pipeline parameters
source('./setup_analysis_pipeline.R')

 tic()
 print(" Collecting results ...")
 source('./script_collect_mcmc_pps_results_multispecies.R')
 print(" done!")
 toc()

