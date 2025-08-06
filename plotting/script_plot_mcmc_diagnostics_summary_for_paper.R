##
# @file script_plot_mcmc_pps_results.R
# @brief plots posterior predictive distributions base on MCMC parameter inference
# using the generalised two-phase recovery model.
#
# @author David J. Warne (david.warne@qut.edu.au)
#           Centre for Data Science
#           School of Mathematical Sciences 
#           Science and Engineering Faculty
#           Queensland University of Technology
#
#
source("setup_analysis_pipeline.R")
library(latex2exp)
library(ggmcmc)
library(ggExtra)
# load filter results with per-capita data 
load(paste(PROC_DATA_DIR,
           sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
           ,".RData",sep=""))
## load MCMC samples
load(paste0(PROC_DATA_DIR,"thinned.conv.chains.ss.v4.RData"))
mcmc.thinned.all.conv.ss <- mcmc.thinned.all.conv


# load spatial data 
gd <- expand.grid(seq(1,10,by=1),seq(5,100,by=5))
gd2 <- expand.grid(seq(0.35,0.75,by=0.1),seq(5,100,by=5))

PMF_T_C0 <- data.frame(T_d=gd[,1],c_0=gd[,2],prob=0)
PMF_alphaD_C0 <- data.frame(alphaD=gd2[,1],c_0=gd2[,2],prob=0)
mcmc.diag.df <- data.frame(Rhat=0,ESS=0)

figsMCMC <- list()
figsMCMC.ss <- list()

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

    mcmc_samples.ss <- mcmc.thinned.all.conv.ss[[k]]

    if (!is.null(mcmc_samples.ss)) { 
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
         print(sprintf("%s site %d %s",reef,site,date))
         diags.ss <- mcmc.diag(mcmc_samples.ss)
         Rhat.ss <- diags.ss$Rhat$mpsrf
         ESS.ss <- min(diags.ss$ESS)
         mcmc.diag.df <- rbind(mcmc.diag.df,c(Rhat.ss,ESS.ss))
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
        
         # look at diagnostics

         samples.ss <- ggs(mcmc_samples.ss)
         levels(samples.ss$Parameter) <- c("alpha","alpha[d]","gamma","T[d]")  
          
        figsMCMC.ss[[k]] <- ggs_traceplot(samples.ss,greek=TRUE) + 
          labs(y="Parameter Value", x="MCMC Iteration",color="Chain")  + 
          scale_colour_brewer(palette = "Dark2")  + 
          theme_bw() + 
          theme(legend.position="top") +
          theme(
            panel.margin.y = unit(-0.5, "lines"),
            strip.background = element_rect(
              color=NA, fill=NA, size=1, linetype="solid"
            ),
            strip.text = element_text(size = 10, face = "italic",vjust=-0.25)
          )
    }
}
mcmc.diag.df <- mcmc.diag.df %>% filter(Rhat != 0, ESS != 0)
mcmc.diag.df.conv <-  mcmc.diag.df %>% filter(Rhat <= 1.1 , ESS >= 200)
mcmc.diag.df.failed <-  mcmc.diag.df %>% filter(Rhat > 1.1 | ESS < 200)

p1 <- ggplot(mcmc.diag.df.conv) + geom_point(aes(x=Rhat,y=ESS)) +labs(x=TeX(r'($\hat{R}$)'),y="ESS") + geom_vline(xintercept = 1.1,linetype="dashed",color="red") + geom_hline(yintercept = 200,linetype="dashed",color="red") + theme_bw()
p1m <- ggMarginal(p1,type="histogram")

p2 <-ggplot(mcmc.diag.df.failed) + geom_point(aes(x=log(Rhat),y=log10(ESS))) +labs(x=TeX(r'($log_{10}(\hat{R})$)'),y=TeX(r'($\log_{10}(ESS)$)')) + geom_vline(xintercept = log10(1.1),linetype="dashed",color="red") + geom_hline(yintercept = log10(200),linetype="dashed",color="red") + theme_bw()
p2m <- ggMarginal(p2,type="histogram")

grid.arrange(grobs = list(p1m,p2m),nrow=1,ncol=2)

 dims <- c(4,5)
# 
inds <- which(!sapply(figsMCMC.ss,is.null))
figsMCMC.ss<-figsMCMC.ss[!sapply(figsMCMC.ss,is.null)]

figsMCMC.ss <- figsMCMC.ss[c(1,3,5,6)]


grid.arrange(grobs = figsMCMC.ss,nrow=1,ncol=4)


