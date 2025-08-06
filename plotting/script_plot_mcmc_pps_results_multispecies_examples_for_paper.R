##
# @file script_plot_mcmc_pps_results.R
# @brief plots posterior predictive distributions base on MCMC parameter inference
# using the generalised two-phase recovery model.
#
# @author David J. Warne (david.warne@qut.edu.au)
#           ARC Centre of Excellence for Mathematical and Statistical Frontiers
#           Centre for Data Science
#           School of Mathematical Sciences 
#           Science and Engineering Faculty
#           Queensland University of Technology
#
#
source("setup_analysis_pipeline.R")
library(latex2exp)
# load filter results with per-capita data 
load(paste(PROC_DATA_DIR,
           sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
           ,".RData",sep=""))
## load MCMC samples
load(paste0(PROC_DATA_DIR,"thinned.conv.chains.ss.v4.RData"))
mcmc.thinned.all.conv.ss <- mcmc.thinned.all.conv
load(paste0(PROC_DATA_DIR,"thinned.conv.chains.v4.RData"))

## load posterior predictive samples
load(paste(PROC_DATA_DIR,"pps.sum.sp.",
           sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
           ,".RData",sep=""))
filt.rec.traj.proc.post.pred.res.ss <- filt.rec.traj.proc.post.pred.res
load(paste(PROC_DATA_DIR,"pps.v4.",
           sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
            ,".RData",sep=""))

# load spatial data 
gd <- expand.grid(seq(1,10,by=1),seq(5,100,by=5))
gd2 <- expand.grid(seq(0.35,0.75,by=0.1),seq(5,100,by=5))

PMF_T_C0 <- data.frame(T_d=gd[,1],c_0=gd[,2],prob=0)
PMF_alphaD_C0 <- data.frame(alphaD=gd2[,1],c_0=gd2[,2],prob=0)

figs <- list()
figsData <- list()
figsT <- list()
figsAD <- list()
figsA <- list()
figsG <- list()

s <- c()
inds <- c()

for (k in 1:730) {
    
    reef  <- filt.rec.traj.proc[[k]]$REEF[1]
    p_code <- filt.rec.traj.proc[[k]]$P_CODE[1]
    depth <- filt.rec.traj.proc[[k]]$DEPTH[1]
    rp_id <- filt.rec.traj.proc[[k]]$RP_ID[1]
    site  <- filt.rec.traj.proc[[k]]$SITE_NO[1]
    date <- lubridate::as_date(filt.rec.traj.proc[[k]]$Date[1])
    s[k] <- filt.rec.traj.proc[[k]]$HC[length(filt.rec.traj.proc[[k]]$HC)]

    mcmc_samples <- mcmc.thinned.all.conv[[k]]
    mcmc_samples.ss <- mcmc.thinned.all.conv.ss[[k]]
    pps <- filt.rec.traj.proc.post.pred.res[[k]]
    pps.ss <- filt.rec.traj.proc.post.pred.res.ss[[k]]
    
    if (!is.null(mcmc_samples) && !is.null(pps) && 
        ((reef == 'Gannett Cay' && site == 2) || #||
         (reef == 'Thetford' && site == 1) ||
          (reef == 'Turner' && site == 1) ||
         (reef == 'Lady Musgrave' && site == 1) 
         )) { 
        inds <- c(inds,k)
        T_data <- c()
        ad_data <- c()
        alpha_data <- c()
        gamma_data <- c()
         T_data_acr <- c()
         alphaD_data_acr <- c()
         alpha_data_acr <- c()
         gamma_data_acr <- c()
         
         # look at diagnostics
         diags.ss <- mcmc.diag(mcmc_samples.ss)
         Rhat.ss <- diags.ss[[1]]$mpsrf
         for (i in 1:length(mcmc_samples.ss)) {
           T_data <- c(T_data,mcmc_samples.ss[[i]][,'Td'])
           ad_data <- c(ad_data,mcmc_samples.ss[[i]][,'alphaD'])
           alpha_data <- c(alpha_data,mcmc_samples.ss[[i]][,'alpha'])
           gamma_data <- c(gamma_data,mcmc_samples.ss[[i]][,'gamma'])
         }
         T_df <- data.frame(T=T_data)
         alphaD_df <- data.frame(alphaD=ad_data)
         alpha_df <- data.frame(alpha=alpha_data)
         gamma_df <- data.frame(gamma=gamma_data)
         nVisits <- ncol(pps.ss)
         pps_stats <- data.frame(colQuantiles(pps.ss[,seq(1,nVisits,by=1)],probs = c(0.005,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.995)))
         pps_stats <- pps_stats %>% mutate(Time = filt.rec.traj.proc[[k]]$T/365.0)
         pps_stats[pps_stats < 0] = 0
         # look at diagnostics
         diags <- mcmc.diag(mcmc_samples)
         Rhat <- diags[[1]]$mpsrf
         print(Rhat)
         for (i in 1:length(mcmc_samples)) {
             T_data_acr <- c(T_data_acr,mcmc_samples[[i]][,'TdA'])
             alphaD_data_acr <- c(alphaD_data_acr,mcmc_samples[[i]][,'alphaDA'])
             alpha_data_acr <- c(alpha_data_acr,mcmc_samples[[i]][,'alphaA'])
             gamma_data_acr <- c(gamma_data_acr,mcmc_samples[[i]][,'gammaA'])
         }
         T_df_acr <- data.frame(T=T_data_acr)
         alphaD_df_acr <- data.frame(alphaD=alphaD_data_acr)
         alpha_df_acr <- data.frame(alpha=alpha_data_acr)
         gamma_df_acr <- data.frame(gamma=gamma_data_acr)
         
         T_data_oth <- c()
         alphaD_data_oth <- c()
         alpha_data_oth <- c()
         gamma_data_oth <- c()
         for (i in 1:length(mcmc_samples)) {
             T_data_oth <- c(T_data_oth,mcmc_samples[[i]][,'TdC'])
             alphaD_data_oth <- c(alphaD_data_oth,mcmc_samples[[i]][,'alphaDC'])
             alpha_data_oth <- c(alpha_data_oth,mcmc_samples[[i]][,'alphaC'])
             gamma_data_oth <- c(gamma_data_oth,mcmc_samples[[i]][,'gammaC'])
         }
         T_df_oth <- data.frame(T=T_data_oth)
         alphaD_df_oth <- data.frame(alphaD=alphaD_data_oth)
         alpha_df_oth <- data.frame(alpha=alpha_data_oth)
         gamma_df_oth <- data.frame(gamma=gamma_data_oth)
         
        nVisits <- ncol(pps)/2
        pps_stats_acr <- data.frame(colQuantiles(pps[,seq(1,nVisits,by=1)],probs = c(0.005,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.995)))
        pps_stats_oth <- data.frame(colQuantiles(pps[,nVisits + seq(1,nVisits,by=1)],probs = c(0.005,0.025,0.05,0.25,0.5,0.75,0.95,0.975,0.995)))
        pps_stats_acr <- pps_stats_acr %>% mutate(Time = filt.rec.traj.proc[[k]]$T/360.0)
        pps_stats_acr[pps_stats_acr < 0] = 0
        pps_stats_oth <- pps_stats_oth %>% mutate(Time = filt.rec.traj.proc[[k]]$T/360.0)
        pps_stats_oth[pps_stats_oth < 0] = 0
        # todo this is a total hack... there must be a better way 
        eval(parse(text = paste0(
        'figs[[k]]  <- ggplot(data=pps_stats_oth) + ggtitle(sprintf("%s site %d %s",reef,site,date)) +
                      geom_ribbon(aes(x=Time,ymin=X0.5.,ymax=X99.5.),alpha=0.25,fill="blue") +
                      geom_ribbon(aes(x=Time,ymin=X2.5.,ymax=X97.5.),alpha=0.25,fill="blue") +
                      geom_ribbon(aes(x=Time,ymin=X5.,ymax=X95.),alpha=0.25,fill="blue") +
                      geom_ribbon(aes(x=Time,ymin=X25.,ymax=X75.),alpha=0.25,fill="blue") +
                      #geom_density(data=T_df_oth, aes(x=T,y=',s[k],'*(..scaled..)),adjust = 2.5, fill = "blue", alpha = 0.3)+
                      geom_ribbon(data=pps_stats_acr,aes(x=Time,ymin=X0.5.,ymax=X99.5.),alpha=0.25,fill="red") +
                      geom_ribbon(data=pps_stats_acr,aes(x=Time,ymin=X2.5.,ymax=X97.5.),alpha=0.25,fill="red") +
                      geom_ribbon(data=pps_stats_acr,aes(x=Time,ymin=X5.,ymax=X95.),alpha=0.25,fill="red") +
                      geom_ribbon(data=pps_stats_acr,aes(x=Time,ymin=X25.,ymax=X75.),alpha=0.25,fill="red") +
                      #geom_density(data=T_df_acr, aes(x=T,y=',s[k],'*(..scaled..)),adjust = 2.5, fill = "red", alpha = 0.3)+
                      geom_ribbon(data=pps_stats,aes(x=Time,ymin=X0.5.,ymax=X99.5.),alpha=0.25,fill="green") +
                      geom_ribbon(data=pps_stats,aes(x=Time,ymin=X2.5.,ymax=X97.5.),alpha=0.25,fill="green") +
                      geom_ribbon(data=pps_stats,aes(x=Time,ymin=X5.,ymax=X95.),alpha=0.25,fill="green") +
                      geom_ribbon(data=pps_stats,aes(x=Time,ymin=X25.,ymax=X75.),alpha=0.25,fill="green") +
                      geom_line(data=filt.rec.traj.proc[[k]],aes(x=T/365,y=HC),linetype = "solid", color="black",size=1) +
                      geom_point(data=filt.rec.traj.proc[[k]],aes(x=T/365,y=HC)) +
                      geom_errorbar(data=filt.rec.traj.proc[[k]],aes(x=T/365,ymin=HC-HC_se,ymax=HC+HC_se)) +
                      geom_line(data=filt.rec.traj.proc[[k]],aes(x=T/365,y=OTH),linetype = "dashed", color="black",size=1) +
                      geom_point(data=filt.rec.traj.proc[[k]],aes(x=T/365,y=OTH)) +
                      geom_errorbar(data=filt.rec.traj.proc[[k]],aes(x=T/365,ymin=OTH-OTH_se,ymax=OTH+OTH_se)) +
                      geom_line(data=filt.rec.traj.proc[[k]],aes(x=T/365,y=ACR),linetype = "dotted", color="black",size=1) +
                      geom_point(data=filt.rec.traj.proc[[k]],aes(x=T/365,y=ACR)) +
                      geom_errorbar(data=filt.rec.traj.proc[[k]],aes(x=T/365,ymin=ACR-ACR_se,ymax=ACR+ACR_se)) +
                      xlim(0,NA) +
                    labs(x="Time (years)",y=" Coral Cover (% Area)") + theme_bw()')))
        eval(parse(text = paste0(
        'figsData[[k]]  <- ggplot(data=pps_stats_oth) + ggtitle(sprintf("%s site %d %s",reef,site,date)) +
                      geom_line(data=filt.rec.traj.proc[[k]],aes(x=T/365,y=HC),linetype = "solid", color="black",size=1) +
                      geom_point(data=filt.rec.traj.proc[[k]],aes(x=T/365,y=HC)) +
                      geom_errorbar(data=filt.rec.traj.proc[[k]],aes(x=T/365,ymin=HC-HC_se,ymax=HC+HC_se)) +
                      geom_line(data=filt.rec.traj.proc[[k]],aes(x=T/365,y=OTH),linetype = "dashed", color="black",size=1) +
                      geom_point(data=filt.rec.traj.proc[[k]],aes(x=T/365,y=OTH)) +
                      geom_errorbar(data=filt.rec.traj.proc[[k]],aes(x=T/365,ymin=OTH-OTH_se,ymax=OTH+OTH_se)) +
                      geom_line(data=filt.rec.traj.proc[[k]],aes(x=T/365,y=ACR),linetype = "dotted", color="black",size=1) +
                      geom_point(data=filt.rec.traj.proc[[k]],aes(x=T/365,y=ACR)) +
                      geom_errorbar(data=filt.rec.traj.proc[[k]],aes(x=T/365,ymin=ACR-ACR_se,ymax=ACR+ACR_se)) +
                      xlim(0,NA) + theme_bw() + scale_x_continuous(breaks=seq(0,11,1)) + scale_y_continuous(breaks=seq(0,100,10)) +
                    labs(x="Time (years)",y=" Coral Cover (% Area)") ')))
        figsT[[k]]  <- ggplot(data=T_df_oth) + 
            geom_density( aes(x=T),adjust = 2.5, fill = "blue", alpha = 0.3) +
            geom_density(data=T_df_acr, aes(x=T),adjust = 2.5, fill = "red", alpha = 0.3) + 
            geom_density(data=T_df, aes(x=T),adjust = 2.5, fill = "green", alpha = 0.3) +
            labs(x=TeX(r'($T_d$)'),y=TeX(r'($p(T_d | Y_{obs})$)')) + theme_bw()
        figsAD[[k]]  <- ggplot(data=alphaD_df_oth) + 
            geom_density( aes(x=alphaD),adjust = 2.5, fill = "blue", alpha = 0.3) +
            geom_density(data=alphaD_df_acr, aes(x=alphaD),adjust = 2.5, fill = "red", alpha = 0.3) + 
          geom_density(data=alphaD_df, aes(x=alphaD),adjust = 2.5, fill = "green", alpha = 0.3) +
            labs(x=TeX(r'($\alpha_d$)'),y=TeX(r'($p(\alpha_d | Y_{obs})$)')) + theme_bw()
        figsA[[k]]  <- ggplot(data=alpha_df_oth) + 
            geom_density( aes(x=alpha),adjust = 2.5, fill = "blue", alpha = 0.3) +
            geom_density(data=alpha_df_acr, aes(x=alpha),adjust = 2.5, fill = "red", alpha = 0.3) + 
            geom_density(data=alpha_df, aes(x=alpha),adjust = 2.5, fill = "green", alpha = 0.3) +
            labs(x=TeX(r'($\alpha$)'),y=TeX(r'($p(\alpha | Y_{obs})$)')) + theme_bw()
        figsG[[k]]  <- ggplot(data=gamma_df_oth) + 
            geom_density( aes(x=gamma),adjust = 2.5, fill = "blue", alpha = 0.3) +
            geom_density(data=gamma_df_acr, aes(x=gamma),adjust = 2.5, fill = "red", alpha = 0.3) + 
            geom_density(data=gamma_df, aes(x=gamma),adjust = 2.5, fill = "green", alpha = 0.3) +
            labs(x=TeX(r'($\gamma$)'),y=TeX(r'($p(\gamma | Y_{obs})$)')) + theme_bw()
    }
}


 dims <- c(4,5)
# 
inds <- which(!sapply(figs,is.null))
figs<-figs[!sapply(figs,is.null)]
figsData<-figsData[!sapply(figsData,is.null)]
figsT<-figsT[!sapply(figsT,is.null)]
figsAD<-figsAD[!sapply(figsAD,is.null)]
figsA<-figsA[!sapply(figsA,is.null)]
figsG<-figsG[!sapply(figsG,is.null)]
figs <- figs[c(1,3,5,6)]
figsData <- figsData[c(1,3,5,6)]
figsT <- figsT[c(1,3,5,6)]
figsAD <- figsAD[c(1,3,5,6)]
figsA <- figsA[c(1,3,5,6)]
figsG <- figsG[c(1,3,5,6)]
inds <- inds[c(1,3,5,6)]

grid.arrange(grobs = figsData,nrow=2,ncol=2)
grid.arrange(grobs = figs,nrow=1,ncol=4)
grid.arrange(grobs = figsT,nrow=1,ncol=4)
grid.arrange(grobs = figsAD,nrow=1,ncol=4)
grid.arrange(grobs = figsA,nrow=1,ncol=4)
grid.arrange(grobs = figsG,nrow=1,ncol=4)


