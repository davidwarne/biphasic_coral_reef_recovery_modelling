#' Model definition of Two-phase recovery for individual multi-species trajectories i
#'
#' @description
#' This file defines the necessary functions to define the \code{model} object for the 
#' single species two-phase model at the level of individual trajectories. This 
#' object is required to use the MCMC functions in and predictive sampling functions
#' from \code{LTMPModellingTools.R}. In effect, this file takes the place of a *.stan 
#' model file, or a JAGS/BUGS model string.
#'
#' @section Warning:
#' As this is cablibrated on individual trajectories, the upper bounds of the \eqn{T_d^A} and \eqn{T_d^C} parameters
#' should be set to the second last observation of the data trajectory (See example pipelines). 
#' As a result, \code{model$upper["TdA"] <- -1; model$upper["TdC"] <- -1}} by default to 
#'force an error if this is not over written.
#' 
#' @note This file (or a copy of it) needs to be modified to change the form of the ODEs
#' the prior definitions and other parameter constraints.
#'
#' @author
#' \itemize{
#'      \item David J. Warne[1,2,3] (\email{david.warne@qut.edu.au})
#' }
#' \enumerate{
#'  \item School of Mathematical Sciences, Faculty of Science, Queensland University of Technology
#'  \item Centre for Data Science, Queensland University of Technology
#'  \item ARC Centre of Excellence for Mathematical and Statistical Frontiers (ACEMS)
#' }
#' 
#' @docType data
#' @name AA_Model_Summary
NULL

# get processed data
load(paste(PROC_DATA_DIR,sprintf(FILTER_REFMT_OUT_DATA_FMT,
                                 max_init,min_final,min_obs)
           ,".RData",sep=""))

# build MCMC sampler configuration
conf <- list(chains = 4,        # number of independent chains to use
            iter = 10000,       # number of sampling interations (or number of interations between diagnostic checks)
            burnin = 10000,     # number of iterations per burnin/warmup step
            CPUs = 4,          # number of CPUs available for parallel chains (optimal CPUs = chains)
            nadapt = 1,        # number of adapatation steps
            initscale = 0.1,
            Rthresh = 1.1,     # stopping criteria threshold for Gelman-Rubin statistic diagnostic check
            ESSthresh = 400,  # stopping criteria threshold for Effective Sample Size diagnostic check
            convcheck = TRUE,
            maxChecks = 300,  # repeat iterations until stopping criteria are satisfied
            maxInits = 10000)


# Perform model calibration to each recovery trajectory using MCMC
tic()
i <- job_id
set.seed(i*1337)
# ensure data is sorted in ascending time order
traj <- filt.rec.traj.proc[[i]] %>% arrange(VISIT_NO)

# build data object
data <- list(nVisits = length(traj$HC[-1]),         # number of visits excluding initial visit
             c0 = traj$HC[1], # cover of initial visit
             t0 = traj$T[1]/365.0,                  # time of initial visit in years
             ts = traj$T[-1]/365.0,                 # time of final visit in years
             K = 100 - traj$AB[length(traj$HC)],    # carrying capacity cover
             C = traj$HC[-1],        # other hard coral cover value time series
             Serr = traj$HC_se) # cover standard error time series

low <- c(0.0,0.0,0.0,0.0)
up <-  c(0.9,1,2.0,traj$T[length(traj$T)]/365)
# build model structure
model <- list(ode_func = general_logistic_twophase,         # RHS for ODE model
              ode_sol = general_logistic_twophase_analytic, 
                loglike = loglike,                          # log likelihood function
                like_sampler = like_sampler,                # simulation of data generation process (for pred. checks)
                logprior = logprior_unif,                   # log prior density
                prior_sampler = prior_sampler_unif,         # prior sampler
                lower = low,                                # lower parameter bounds 
                upper = up,                                 # upper parameter bounds
                hyp = c(low[1],up[1],low[2],up[2],low[3],up[3],low[4],up[4]),   # prior hyper-parameters
                varnames = vlab)                            # parameter labels
fprintf("MCMC sampling for recovery trajectory %d\n",i)
tic()
samples <- adaptMCMC_fit_ode_model(data,model,conf)
toc()

filt.rec.traj.proc.mcmc.res <- samples
fprintf("Posterior predictive sampling for recovery trajectory %d\n",i)
tic()
filt.rec.traj.proc.post.pred.res <- posterior_predictive_samples(data,model,samples)
toc()
#
# # save results
save(filt.rec.traj.proc.mcmc.res,file=paste(PROC_DATA_DIR,"mcmc.sum.sp.",
                   sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
                                        ,".jobid",num2str(i),".RData",sep=""))

save(filt.rec.traj.proc.post.pred.res,file=paste(PROC_DATA_DIR,"pps.sum.sp.",
                   sprintf(FILTER_REFMT_OUT_DATA_FMT,max_init,min_final,min_obs)
                                        ,".jobid",num2str(i),".RData",sep=""))
