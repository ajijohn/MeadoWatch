# Phenology functions
# @Authors Janneke 

#Null model - assumes the probability of flowering does not vary with time
nullfit <- function (param){ 
  meanp  <- param[1]
  pred   <- rep(meanp, times=length(phenophase))
  llik     <- dbinom(phenophase,1,pred, log=TRUE)
  return(-sum(llik))
}

#Alternative model - assumes flowering varies with DOY according to a guassian curve
#can fit this to individual plot / year / species data, data for a species in all plots, etc
#curve fitting function for mle
curvefit <- function (param){ 
  peakp  <- param[1]
  rangep <- param[2]
  maxp   <- param[3]
  #pred   <- maxp*(1/(rangep*sqrt(2*pi)))*exp(-0.5*((days-peakp)/rangep)^2)
  pred <- inv.logit(rangep * (days - peakp)^2 + maxp)
  llik     <- dbinom(phenophase,1,pred, log=TRUE)
  return(-sum(llik))
}

#Fits model per species 
#uses snow disappearance as predictor of peak flowering, plot/trail specific range and
curvefit_snowmod <- function (param){ #curve fitting function for mle
  intpeak  <- param[1]
  slopepeak <- param[2]
  rangep_all <- param[3:(2+ntrlyr)] #different range per trail, year
  maxp_all   <- param[(3+ntrlyr):(2+2*ntrlyr)] #different max per trail, year
  peakp <- intpeak + slopepeak*SDD
  rangep <- rangep_all[trlyr]
  maxp <- maxp_all[trlyr]
  #pred   <- maxp*(1/(rangep*sqrt(2*pi)))*exp(-0.5*((days-peakp)/rangep)^2)
  pred <- inv.logit(rangep * (days - peakp)^2 + maxp)
  llik     <- dbinom(phenophase,1,pred, log=TRUE)
  return(-sum(llik))
}


#Predicts probability of observing phenological phase, given a vector of days and parameters (used for plotting)
predphen <- function (xx, param){ 
  days    <- xx
  peakp  <- param[1]
  rangep <- param[2]
  maxp   <- param[3]
  #pred   <- maxp*(1/(rangep*sqrt(2*pi)))*exp(-0.5*((days-peakp)/rangep)^2)
  pred <- inv.logit(rangep * (days - peakp)^2 + maxp)
  return(pred)
}