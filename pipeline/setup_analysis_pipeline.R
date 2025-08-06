##
# @file setup_analysis_pipeline.R
# @brief complete initial pipeline for reproducibility.
#
# @author David J. Warne (david.warne@qut.edu.au)
#           ARC Centre of Excellence for Mathematical and Statistical Frontiers
#           School of Mathematical Sciences 
#           Faculty of Science
#           Queensland University of Technology
# 
#

library(tidyverse)
library(ggmap)
library(gridExtra)
library(functional)
library(tictoc)
library(pracma)
library(gap)
library(stats4)
library(scales)
library(ggExtra)
library(cowplot)
library(matrixStats)

library(mcmc)
library(adaptMCMC)
library(postpack)
library(coda)
library(deSolve)
library(parallel)
library(doRNG)
library(doFuture)
options(future.globals.maxSize = 8000 * 1024^2)
source("../LTMPTools/LTMPDataTools.R")
source("../LTMPTools/LTMPModellingTools.R")

# input dir and data
SCRIPT_DIR                <- "../"
MODEL_DIR                 <- "../models/"
DATA_DIR                  <- "../data/primary/"
DIST_DAT_FILE             <- "disturbance.RData"
SAMPLE_DAT_FILE           <- "samples.RData"
VISIT_DAT_FILE            <- "ker.code.visit.RData"
SPAT_DAT_FILE             <- "spatial.dat.RData"
GRP_LEV_TRANS_DAT_FILE    <- "groups.transect.RData"
GRP_ACR_LEV_TRANS_DAT_FILE <- "groups.acr.transect.RData"
BNTHS_LEV_TRANS_DAT_FILE  <- "benthos.transect.RData"
FAMILY_LEV_TRANS_DAT_FILE  <- "family.cover.transect.RData"


# output dir and data
PROC_DATA_DIR           <- "./temp_processed/"
REC_TRAJ_DAT_FILE       <- "time.series.site.transect.RData"
TIME_SERIES_DAT_FILE    <- "recovery.trajectories.site.transect.RData"

FILTER_OUT_DATA_FMT         <- "rec.traj.trans.init%f.final%f.obs%d"
FILTER_REFMT_OUT_DATA_FMT   <- "rec.traj.trans.proc.init%f.final%f.obs%d"

# for any automatic figure outputs
FIG_DIR <- "../data/figures/"

# filter critera 
max_init    <- 100.0
min_final   <- 0.0
min_obs     <- 3.0
