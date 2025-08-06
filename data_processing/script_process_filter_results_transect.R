##
# @file script_process_filter_results.R
# @brief script to transform recovery trajectories into a more suitable form for
# individual analysis.
#
# @author David J. Warne (david.warne@qut.edu.au)
#           ARC Centre of Excellence for Mathematical and Statistical Frontiers
#           School of Mathematical Sciences 
#           Science and Engineering Faculty
#           Queensland University of Technology
#
#

# load filter results 
load(paste(PROC_DATA_DIR,sprintf(FILTER_OUT_DATA_FMT,max_init,min_final,min_obs)
           ,".RData",sep=""))

filt.rec.traj.proc <- list()

j <- 1
# reformat all trajectories and store in a list
for (rp_id in unique(filt.rec.traj$RP_ID)) {
   cover.dat <- filt.rec.traj %>% filter(RP_ID == rp_id)
   filt.rec.traj.proc[[j]] <- reformat_recovery_trajectories(cover.dat)
   j <- j + 1
}

# save the cleaned-up trajectory data
save(filt.rec.traj.proc,file=paste(PROC_DATA_DIR,
                  sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
                  ,".RData",sep=""))
