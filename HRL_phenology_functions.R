# Phenology functions
# Used for MeadoWatch / WTA markdown 
# Janneke Hille Ris Lambers
# December 15, 2020

# Null model - assumes the probability of flowering is constant
# this is a binomial likelihood model that estimates a constant probability
# of flowering (one parameter) and the negative log likelihood, 
# given yes / no observations. This is a useful null model 
# to test the phenological stage of interest is unimodal.

nullfit <- function (param){ 
  meanp  <- param #mean probability of observing phenological phase
  pred   <- rep(meanp, times=length(phenophase))
  llik     <- dbinom(phenophase,1,pred, log=TRUE) #this is the likelihood
  return(-sum(llik)) #this is the negative log likelihood, which is minimized
}


# Fits a phenological curve per plot. A binomial likelihood model that 
# estimates the probability of flowering as a function of DOY; as described
# by three parameters - range, maximum, and peak. The unimodal curve
# is fit with a logit transformation (similar to Theobald et al 2017, 
# Sethi et al 2020; albeit with parameters that do NOT vary with climate). 
# Also requires yes / no observations. This function can be fit to 
# yes / no observations of phenology (phenophase) observed by DOY

curvefit_perplot <- function (param){ #curve fitting function for mle
  peakp  <- param[1] #in days, the time of peak phenophase
  rangep <- param[2] #parameter describing width of phenological curve
  maxp   <- param[3] #parameter describing height of phenological curve
  pred <- inv.logit(rangep * (days - peakp)^2 + maxp)
  llik     <- dbinom(phenophase,1,pred, log=TRUE)
  return(-sum(llik))
}


# Phenological mle model per species, fit to yes / no data from multiple plots
# Also a binomial likelihood model that assuming that flowering probability
# for a species in a plot is a function of several parameters. The probability 
# of flowering is assumed to be unimodal with peak flowering in a plot 
# depending on snowmelt date at that plot (a simple linear relationship
# described by an intercept and slope parameter). The parameters describing
# the range and maximum are also fit per species

curvefit_species <- function (param){ #curve fitting function for mle
  intpeak  <- param[1]
  slopepeak <- param[2]
  rangep <- param[3]
  maxp <- param[4]
  peakp <- intpeak + slopepeak*SDD
  pred <- inv.logit(rangep * (days - peakp)^2 + maxp)
  llik     <- dbinom(phenophase,1,pred, log=TRUE)
  return(-sum(llik))
}


# Estimate peak flowering as a function of peak, duration and max parameters.
# This function takes as input DOY (days since January 1) and 3 parameters 
# describing peak, range and duration of flowering; and returns the 
# probability of flowering (primarily used for plotting). 

predflower <- function (xx, param){ 
  days    <- xx
  peakp  <- param[1]
  rangep <- param[2]
  maxp   <- param[3]
  pred <- inv.logit(rangep * (days - peakp)^2 + maxp)
  return(pred)
}

# Estimate the probability of observing a flower within a MeadoWatch plot.
# This function takes as input DOY (days since January 1) and 4 parameters 
# describing the earliest peak flowering observed along the trail, the latest
# peak flowering observed along the trail, and parameters describing the duration
# and maximum probability of flowering. 

predflowertrail <- function (xx, param){ 
  days    <- xx
  peakearly  <- param[1]
  peaklate  <- param[2]
  range <- param[3]
  max   <- param[4]
  predmax <- inv.logit(max)
  ptrailflower <- rep(predmax, time=length(days))
  predearly <- inv.logit(range * (days - peakearly)^2 + max)
  predlate <- inv.logit(range * (days - peaklate)^2 + max)
  ptrailflower[days[]<peakearly] <- predearly[days[]<peakearly]
  ptrailflower[days[]>peaklate] <- predlate[days[]>peaklate]
  ptrailflower[ptrailflower<0.005] <- 0
  return(ptrailflower)
}


#########FOR POSTERITY - no longer in use

#Fits model per species 
#peak flowering varies with snowmelt, and plot/trail specific range and max
#Should be fit per species (but all data for each speices)
# A binomial likelihood model that estimates the probability of flowering
# for a species across all plots, years and transects as a function of 
# several parameters. The probability of flowering is assumed to be unimodal
# with peak flowering depending on snowmelt date at that plot (a simple linear
# relationship described by an intercept and slope parameter). Parameters 
# describing the range and maximum flowering are allowed to vary by year and
# and trail. This function is fit to all data from a species.
curvefit_allplot <- function (param){ #curve fitting function for mle
  intpeak  <- param[1]
  slopepeak <- param[2]
  rangep_all <- param[3:(2+ntrlyr)] #different range per trail, year
  maxp_all   <- param[(3+ntrlyr):(2+2*ntrlyr)] #different max per trail, year
  peakp <- intpeak + slopepeak*SDD
  rangep <- rangep_all[trlyr]
  maxp <- maxp_all[trlyr]
  pred <- inv.logit(rangep * (days - peakp)^2 + maxp)
  llik     <- dbinom(phenophase,1,pred, log=TRUE)
  return(-sum(llik))
}
