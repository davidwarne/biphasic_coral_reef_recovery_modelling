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
library(ggpmisc)
# load filter results with per-capita data 
load(paste(PROC_DATA_DIR,
           sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
           ,".RData",sep=""))
## load MCMC samples
load(paste0(PROC_DATA_DIR,"thinned.conv.chains.ss.v4.RData"))
## load posterior predictive samples
load(paste(PROC_DATA_DIR,"pps.sum.sp.",
           sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
            ,".RData",sep=""))
# load spatial data 
gd <- expand.grid(seq(1,10,by=1),seq(5,100,by=5))
gd2 <- expand.grid(seq(0.35,0.75,by=0.1),seq(5,100,by=5))

PMF_T_C0 <- data.frame(T_d=gd[,1],c_0=gd[,2],prob=0)
PMF_alphaD_C0 <- data.frame(alphaD=gd2[,1],c_0=gd2[,2],prob=0)


figs <- list()
figsT <- list()
figsAD <- list()
figsA <- list()
figsG <- list()
pps.dat.df <- data.frame(COVER=0,quantile=0,CrI=0)
pred.vs.obs.df <- data.frame(PRED=0,OBS=0)
s <- c()
inds <- c()
kk <- 1
failed <- c()
# test thresholds
Td_crit = 2.0
alphaD_crit = 0.75
for (k in 1:737) {
    reef  <- filt.rec.traj.proc[[k]]$REEF[1]
    p_code <- filt.rec.traj.proc[[k]]$P_CODE[1]
    depth <- filt.rec.traj.proc[[k]]$DEPTH[1]
    rp_id <- filt.rec.traj.proc[[k]]$RP_ID[1]
    site  <- filt.rec.traj.proc[[k]]$SITE_NO[1]
    dyear <- year(filt.rec.traj.proc[[k]]$Date[1])
    Tmax <- max(filt.rec.traj.proc[[k]]$T)/365
    s[k] <- filt.rec.traj.proc[[k]]$HC[length(filt.rec.traj.proc[[k]]$HC)]
    nobs <- nrow(filt.rec.traj.proc[[k]])
    if (k > length(mcmc.thinned.all.conv) ){
        mcmc_samples <- NULL
        pps <- NULL
    }else{
        mcmc_samples <- mcmc.thinned.all.conv[[k]]
        pps <- filt.rec.traj.proc.post.pred.res[[k]]
    }
    if (!is.null(mcmc_samples) && !is.null(pps)) {
      
      inds <- c(inds,k)
        T_data <- c()
        ad_data <- c()
        alpha_data <- c()
        gamma_data <- c()
        # look at diagnostics
        diags <- mcmc.diag(mcmc_samples)
        Rhat <- diags$Rhat$mpsrf
        ESS <- min(diags$ESS)
        if (Rhat <= 1.1 & ESS >= 200){
        
        
          nVisits <- ncol(pps)
          for (i in 1:nVisits){
            dat <- filt.rec.traj.proc[[k]]$HC[i]
            quant <- mean(dat <= pps[,i])
            CrI <- 0
            if (quant < 0.5) {
              CrI <- 1-2*quant
            } else {
              CrI <- 1-2*(1-quant)
            }
            pps.dat.df <- rbind(pps.dat.df,c(dat,quant,CrI))
          }
          
          for (i in 1:nVisits){
            obs <- filt.rec.traj.proc[[k]]$HC[i]*ones(2,1)
            pred <- sample(pps[,i],2,replace = TRUE)
            pred[pred < 0] <-0
            pred[pred > 100] <- 100
            pred.vs.obs.df <- rbind(pred.vs.obs.df,data.frame(PRED=pred,OBS=obs))
            
          }
        }
        
     } 
}


beta_exp <- seq(0.01,0.99,by=0.01)
beta_obs <- 0*beta_exp
beta_obs_var <- 0*beta_exp
beta_obs_std <- 0*beta_exp
for (i in 1:length(beta_exp)){ 
  beta_obs[i] <- sum(pps.dat.df$CrI <= beta_exp[i])/nrow(pps.dat.df)
  beta_obs_var[i] <- beta_obs[i]*(1 - beta_obs[i])/nrow(pps.dat.df)
  beta_obs_std[i] <- sqrt(beta_obs_var[i])
}
beta.df <- data.frame(beta_exp=100*beta_exp,beta_obs=100*beta_obs,beta_std=100*beta_obs_std)

q_exp <- seq(0.0,1,by=0.001)
q_obs <- 0*q_exp
CrI <- 0*q_exp
for (i in 1:length(q_exp)){ 
  q_obs[i] <- sum(pps.dat.df$quantile <= q_exp[i])/nrow(pps.dat.df) 
}
qq.df <- data.frame(q_exp=q_exp,q_obs=q_obs)
p1 <- ggplot(data=beta.df) + 
  geom_errorbar(aes(x=beta_exp,ymin=beta_obs-1.96*beta_std,ymax=beta_obs+1.96*beta_std)) +
  geom_line(aes(x=beta_exp,y=beta_obs)) + 
  geom_abline(slope = 1,intercept = 0,color="red") + 
  labs(x=TeX(r'($\beta$)'),y=TeX(r'(Proportion of data within ${\beta\%}CrI$)')) + 
  scale_y_continuous(limits=c(0,100),expand = c(0,0)) + 
  scale_x_continuous(limits=c(0,100),expand = c(0,0))+
  theme_bw()
p2 <- ggplot(data=pps.dat.df,aes(x=100*CrI,y=COVER)) + 
  geom_density2d_filled(show.legend = FALSE) +

  labs(x=TeX(r'(Smallest Contaning CrI)'),y="Coral Cover (% Area)") + 
  scale_y_continuous(limits=c(0,100),expand = c(0,0)) + 
  scale_x_continuous(limits=c(0,100),expand = c(0,0))+
  theme_bw()

p3 <- ggplot(data=pred.vs.obs.df) +
  geom_point(aes(x=OBS,y=PRED),alpha = 0.1) +
  geom_abline(slope = 1,intercept = 0,color="red") + 
  scale_y_continuous(limits=c(0,100),expand = c(0,0)) + 
  scale_x_continuous(limits=c(0,100),expand = c(0,0))+
  labs(x=TeX(r'(Observed Coral Cover (% Area))'),y=TeX(r'(Predicted Coral Cover (% Area))')) +
  theme_bw()

grid.arrange(grobs =list(p3,p1,p2),nrow=1,ncol=3)

inds <- which(!sapply(figs,is.null))
figs<-figs[!sapply(figs,is.null)]
figsT<-figsT[!sapply(figsT,is.null)]
figsAD<-figsAD[!sapply(figsAD,is.null)]
figsA<-figsA[!sapply(figsA,is.null)]
figsG<-figsG[!sapply(figsG,is.null)]
figs <- figs[c(1,3,5,6)]
figsT <- figsT[c(1,3,5,6)]
figsAD <- figsAD[c(1,3,5,6)]
figsA <- figsA[c(1,3,5,6)]
figsG <- figsG[c(1,3,5,6)]
inds <- inds[c(1,3,5,6)]


grid.arrange(grobs = figs,nrow=1,ncol=4)
grid.arrange(grobs = figsT,nrow=1,ncol=4)
grid.arrange(grobs = figsAD,nrow=1,ncol=4)
grid.arrange(grobs = figsA,nrow=1,ncol=4)
grid.arrange(grobs = figsG,nrow=1,ncol=4)

n <- ceil(length(figs)/25)
for (i in 1:n){
    print(i)
    pdf(file=sprintf('pps_marginals_ss_class_%d.pdf',i),width = 25, height = 25)
    grid.arrange(grobs = figs[seq((i-1)*25 +1,i*25,by=1)],nrow=5,ncol=5)
    dev.off()
}


