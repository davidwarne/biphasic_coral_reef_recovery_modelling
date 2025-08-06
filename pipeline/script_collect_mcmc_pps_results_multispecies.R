##
# @file script_main_init_analysis_pipeline.R
# @brief complete initial pipeline for reproducibility.
#
# @author David J. Warne (david.warne@qut.edu.au)
#           Faculty of Science
#           Queensland University of Technology
# 
#
source('./setup_analysis_pipeline.R')
# create structures for mcmc runs and posterior predictive samples
mcmc.res.all <- list();
mcmcFailed <- c()
pps.res.all <- list();
ppsFailed <- c()
# loop over all (return list of those for which not output produced)
for (i in 1:737) {
    ppsFile <- paste0(PROC_DATA_DIR,"pps.sum.sp.",
        sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs),
        ".jobid",num2str(i),".RData")
    if (file.exists(ppsFile)) {
        load(ppsFile)
        pps.res.all[[i]] <- filt.rec.traj.proc.post.pred.res
    } else {
        ppsFailed <- c(ppsFailed,i)
    }
}
filt.rec.traj.proc.post.pred.res <- pps.res.all;

save(filt.rec.traj.proc.post.pred.res,file=paste(PROC_DATA_DIR,"pps.sum.sp.",
                   sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
                                        ,".RData",sep=""))
