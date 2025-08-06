#' Model definition of Two-phase recovery for individual multi-species trajectories 
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
#' As a result, \code{model$upper["TdA"] <- -1; model$upper["TdC"] <- -1} by default to 
#' force an error if this is not over written.
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

#' Parameter labels
#'
#' @description Returns parameter vector with labels 
#'
#' @note The function relies upon a prespecified list of strings call \code{param_names}. 
#' 
#' @param theta parameter vector
#'
#' @return parameter vector with labels
#'
vlab <- function(theta) {
    labs <- c(1:length(theta))
    k <- 1
    for (L in param_names) {
        labs[k] <- L
        k <- k + 1
    }
    names(theta) <- labs
    return(theta)
}

#' Log density of the prior
#'
#' @description log density of the prior \eqn{p(\theta | \lambda)} where \epn{\lambda} is a vector of hyper-parameters
#'
#' @param theta Parameter vector
#' @param lambda Hyper-parameter vector
#'
logprior <- function(theta,lambda) {
    d <- sum(c(dnorm(theta[1], mean = lambda[1], sd = lambda[2], log = TRUE),
               dnorm(theta[2], mean = lambda[3], sd = lambda[4], log = TRUE),
               dnorm(theta[3], mean  = lambda[5], sd = lambda[6], log = TRUE),
               dnorm(theta[4], mean = lambda[7], sd = lambda[8], log = TRUE),
               dnorm(theta[5], mean = lambda[1], sd = lambda[2], log = TRUE),
               dnorm(theta[6], mean = lambda[3], sd = lambda[4], log = TRUE),
               dnorm(theta[7], mean  = lambda[5], sd = lambda[6], log = TRUE),
               dnorm(theta[8], mean = lambda[7], sd = lambda[8], log = TRUE)))
    if(is.nan(d) || is.na(d)){
        return(-Inf)
    } else {
        return(d)
    }
}
logprior_unif <- function(theta,lambda) {
    d <- sum(c(dunif(theta[1], min = lambda[1], max = lambda[2], log = TRUE),
               dunif(theta[2], min = lambda[3], max = lambda[4], log = TRUE),
               dunif(theta[3], min  = lambda[5], max = lambda[6], log = TRUE),
               dunif(theta[4], min = lambda[7], max = lambda[8], log = TRUE),
               dunif(theta[5], min = lambda[1], max = lambda[2], log = TRUE),
               dunif(theta[6], min = lambda[3], max = lambda[4], log = TRUE),
               dunif(theta[7], min  = lambda[5], max = lambda[6], log = TRUE),
               dunif(theta[8], min = lambda[7], max = lambda[8], log = TRUE)))
    if(is.nan(d) || is.na(d)){
        return(-Inf)
    } else {
        return(d)
    }
}

#' Prior sampler
#'
#' @description Prior sampler for chain initialisation by generateion \eqn{\theta \sim p(.|\lambda)}. 
#' @note In general this function need not actually sample the prionr, but it at least must be 
#' consistent with log-prior is terms of support. 
#' To this end, the function can be used to specify many initialisation strategies
#'
#' @param lambda hyper-parameter vector
prior_sampler <- function(lambda) {
    theta <- c(rnorm(1, mean = lambda[1], sd = lambda[2]),
               rnorm(1, mean = lambda[3], sd = lambda[4]),
               rnorm(1, mean  = lambda[5], sd = lambda[6]),
               rnorm(1, mean = lambda[7], sd = lambda[8]),
               rnorm(1, mean = lambda[1], sd = lambda[2]),
               rnorm(1, mean = lambda[3], sd = lambda[4]),
               rnorm(1, mean  = lambda[5], sd = lambda[6]),
               rnorm(1, mean = lambda[7], sd = lambda[8])) 
    return(theta)
}

prior_sampler_unif <- function(lambda) {
    theta <- c(runif(1, min = lambda[1], max = lambda[2]),
               runif(1, min = lambda[3], max = lambda[4]),
               runif(1, min  = lambda[5], max = lambda[6]),
               runif(1, min = lambda[7], max = lambda[8]),
               runif(1, min = lambda[1], max = lambda[2]),
               runif(1, min = lambda[3], max = lambda[4]),
               runif(1, min  = lambda[5], max = lambda[6]),
               runif(1, min = lambda[7], max = lambda[8])) 
    return(theta)
}


#' Multispecies generalised two-phase model
#'
#' @description Implements right-hand-side (RHS) of the ODE model for coral recovery
#' @note This function is in a specific format for the \code{deSolve} package.
general_logistic_twophase <- function(t, X, theta) {
    with(as.list(c(X,theta)), {
        if (t < TdA) {
            dAdt <- alphaDA*(alphaA/gammaA)*A*(1.0 - ((A + C)/K)^(gammaA))
            #dAdt <- alphaDA*(alphaA/gammaA)*A
        } else {
            dAdt <-         (alphaA/gammaA)*A*(1.0 - ((A + C)/K)^(gammaA))
        }
        if (t < TdC) {
            dCdt <- alphaDC*(alphaC/gammaC)*C*(1.0 - ((A + C)/K)^(gammaC))
            #dCdt <- alphaDC*(alphaC/gammaC)*C
        } else {
            dCdt <-         (alphaC/gammaC)*C*(1.0 - ((A + C)/K)^(gammaC))
        }
        list(c(dAdt,dCdt)) # A = Achropodiae, C = HC (hard corals) - A,
    })
}

#' Log likelihood
#'
#' @description Log likelihood for the data given parameters theta under the model 
#'
#' @param data Data list 
#' @param model Model object 
#' @param theta Parameter vector
#'
loglike <- function(data,model,theta) {
    #extract initial condition (a vector)
    X0 = data$c0
    # labelled parameter vector (with carrying capacity, K, appended)
    p <- c(model$varnames(theta), K = data$K)
    #solve forwards problem for given theta
    X_mu <- ode(y = X0, times = c(data$t0,data$ts), 
                func = model$ode_func, parms = p, method = "ode45")
    d <- 0
    # observation error model
    # set 0 variances to really small non-zeros to avoid nans and bad fits
    sdA <- data$Serr[-1,'A'] 
    sdC <- data$Serr[-1,'C']
    sdA[sdA <= 0] <- 0.1
    sdC[sdC <= 0] <- 0.1
    d <- d + sum(dnorm(data$A,mean = X_mu[-1,'A'],sd = sdA, log = TRUE))    
    d <- d + sum(dnorm(data$C,mean = X_mu[-1,'C'],sd = sdC, log = TRUE))
    
    if (is.nan(d) || is.na(d)) {
        return(-Inf)
    } else {
        return(d)
    }
}

#' Likelihood sampler
#' 
#' @description Generates synthetic data given a fix parameter value \eqn{\theta}. 
#' This is not used by the MCMC functions, but the posterior and prior predictive samplers.
#'
#' @param data Data list
#' @param model Model Object
#' @param theta Parameter vector
like_sampler <- function(data,model,theta) {
    #print(data)
    #extract initial condition
    X0 = data$c0
    # labelled parameter vector (with carrying capacity, K, appended)
    p <- c(theta, K = data$K)
    #solve forwards problem for given theta
    X_mu <- ode(y = X0, times = c(data$t0,data$ts), 
                func = model$ode_func, parms = p, method = "ode45")
    sdA <- data$Serr[,'A'] 
    sdC <- data$Serr[,'C']
    sdA[sdA <= 0] <- 0.1
    sdC[sdC <= 0] <- 0.1
    # observation error model
    return(c(rnorm(length(X_mu[,'A']), mean = X_mu[,'A'], sd = sdA),
             rnorm(length(X_mu[,'C']), mean = X_mu[,'C'], sd = sdC)))
}
# ROXYGEN_STOP
    ## define parameter names
    param_names <- c("alphaDA","alphaA", "gammaA","TdA","alphaDC","alphaC", "gammaC","TdC")
    # build model structure
    model <- list(ode_func = general_logistic_twophase,       # RHS for ODE model
                    loglike = loglike,                          # log likelihood function
                    like_sampler = like_sampler,                # simulation of data generation process (for pred. checks)
                    logprior = logprior,                        # log prior density
                    prior_sampler = prior_sampler,              # prior sampler
                    lower = c(0.25,0,0,0,0.25,0,0,0),                      # lower parameter bounds 
                    upper = c(0.75,Inf,Inf,-1,0.75,Inf,Inf,-1),  # upper parameter bounds
                    hyp = c(0.5,0.5,0.3,0.3,0.0,0.5,0.0,4.0),   # prior hyper-parameters
                    varnames = vlab)                            # parameter labels
  
# ROXYGEN_START
