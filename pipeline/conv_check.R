##
# @file conv_check.R
# @brief Collecting info on convergence diangostics to schedule reruns


convtable <- data.frame(JobID=numeric(),
                    Rhat=numeric(),
                    ESS=numeric(),
                    conv=logical())



mcmc.thinned.all.conv <- list()

job_vec=c()
rhat_vec=c()
ESS_vec=c()
conv_vec=c()

toc_tot <- 0
for (k in 1:737) {
    tic()
    # load MCMC trace for this trajectory
    mcmcFilek <- paste(PROC_DATA_DIR,"mcmc.sum.pcHC.",
           sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs),
           ".jobid",k,".000.RData",sep="") 
    print(mcmcFilek)
    job_vec <- c(job_vec,k)
    if (file.exists(mcmcFilek)){
        load(mcmcFilek)
        
        diags <- mcmc.diag(filt.rec.traj.proc.mcmc.res)
        print(Rhatk <- diags[[1]]$mpsrf)
        Rhat <- diags[[1]]$mpsrf
        ESS <- min(diags[[2]])
        rhat_vec <- c(rhat_vec,Rhat)
        ESS_vec <- c(ESS_vec,ESS)
        conv_vec <- c(conv_vec,Rhat <= 1.1 && ESS >= 400)
        print(c(k,Rhat,ESS,Rhat <= 1.1 && ESS >= 400))
        
        if (Rhat <= 1.1 && ESS >= 400) {
            mcmc.pd <-post_dim(filt.rec.traj.proc.mcmc.res)
            print(mcmc.pd)
            filt.rec.traj.proc.mcmc.res.thin <- post_thin(filt.rec.traj.proc.mcmc.res, keep_iters = ESS)
            print(post_dim(filt.rec.traj.proc.mcmc.res.thin))
            mcmc.diag(filt.rec.traj.proc.mcmc.res.thin)
            mcmc.thinned.all.conv[[k]] <- filt.rec.traj.proc.mcmc.res.thin
        }
        
        # clean up to stop memory leakage
        rm(file.rec.traj.proc.mcmc.res)
        gc() # gc() has side effect of invoking garbage collector 
    } else {
        rhat_vec <- c(rhat_vec,-1)
        ESS_vec <- c(ESS_vec,-1)
        conv_vec <- c(conv_vec,0)
    }
    toc_k <- toc()
    toc_tot <- toc_tot + toc_k
    print(c((toc_tot/k)*(737-k)/3600, 100*k/737))
}

convtable <- data.frame(JobID=job_vec,
                        Rhat=rhat_vec,
                        ESS=ESS_vec,
                        conv=conv_vec)
save(convtable,file="conv.table.v4.RData")
save(mcmc.thinned.all.conv,file="thinned.conv.chains.v4.RData")
