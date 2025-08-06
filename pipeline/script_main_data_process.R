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

# step 1: extract all recovery trajectories (only needed once)
tic()
print(" Extracting recovery trajectories...")
source(paste0(SCRIPT_DIR,'data_processing/script_extract_recovery_trajectories_transect_multispecies.R'))
print(" done!")
toc()

# step 2a: filter on selected criteria
tic()
print(" Applying filter ...")
source(paste0(SCRIPT_DIR,'data_processing/script_filter_recovery_trajectories_transect_multispecies.R'))
print(" done!")
toc()

# step 2b: reformat filter results
tic()
print(" Reformatting filter results ...")
source(paste0(SCRIPT_DIR,'data_processing/script_process_filter_results_transect.R'))
print(" done!")
toc()

