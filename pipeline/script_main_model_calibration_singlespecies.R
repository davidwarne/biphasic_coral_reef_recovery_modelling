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

# build model
source(paste0(MODEL_DIR,'DefineTwoPhaseGeneralModelSingleSpecies.R'))

tic()
print(" Calibrate model and Posterior predictive checks ...")
fprintf(" Job id = %d\n", job_id)
source('./script_apply_mcmc_general_model_singlespecies.R')
print(" done!")
toc()

