## Author: Meera Lee Sethi
## Last worked on: 12-27-2019
## Script for setting up the Bayesian hierarchical model used in "Early snowmelt and warmer, drier summers shrink post-flowering transition times in subalpine wildflowers"

## Sets up workspace 
library(runjags)

##### Reads in phenodata:
train_data <- read.csv("./final_train_data.csv", header=TRUE)

### Writes model/s
 
logit <- function(x){log(x/(1-x))} #Converts from a scale of 0-1 to a scale of -infinity to +infinity

##### 
#####JAGS Model 
 
cat("
    model{
    # priors for community mean shape of phenological curves
    height.mu ~ dnorm(2,0.01)
    height.sigma ~ dunif(0.001,10)
    height.tau <- pow(height.sigma,-2)
    peakflower.sigma ~ dunif(0,10)
    peakflower.tau <- pow(peakflower.sigma,-2)
    flowertoseed.sigma ~ dunif(0,2)
    flowertoseed.tau <- pow(flowertoseed.sigma,-2)
    width.sigma ~ dunif(0,1)
    width.tau <- pow(width.sigma,-2)
    
    
    
    # priors for community mean environmental effects
    peakflower.SDDDOY.mu ~ dnorm(1,0.1)T(0,2)
    peakflower.SDDDOY.sigma ~ dunif(0.001,10)
    peakflower.SDDDOY.tau <- pow(peakflower.SDDDOY.sigma,-2)
    peakflower.gddearly.mu ~ dnorm(0,0.1)
    peakflower.gddearly.sigma ~ dunif(0.001,10)
    peakflower.gddearly.tau <- pow(peakflower.gddearly.sigma,-2)
    peakflower.gddlate.mu ~ dnorm(0,0.1)
    peakflower.gddlate.sigma ~ dunif(0.001,10)
    peakflower.gddlate.tau <- pow(peakflower.gddlate.sigma,-2)
    peakflower.moisture.mu ~ dnorm(0,0.1)
    peakflower.moisture.sigma ~ dunif(0.001,10)
    peakflower.moisture.tau <- pow(peakflower.moisture.sigma,-2)
    flowertoseed.SDDDOY.mu ~ dnorm(0,0.1)
    flowertoseed.SDDDOY.sigma ~ dunif(0.001,10)
    flowertoseed.SDDDOY.tau <- pow(flowertoseed.SDDDOY.sigma,-2)
    flowertoseed.gddearly.mu ~ dnorm(0,0.1)
    flowertoseed.gddearly.sigma ~ dunif(0.001,10)
    flowertoseed.gddearly.tau <- pow(flowertoseed.gddearly.sigma,-2)
    flowertoseed.gddlate.mu ~ dnorm(0,0.1)
    flowertoseed.gddlate.sigma ~ dunif(0.001,10)
    flowertoseed.gddlate.tau <- pow(flowertoseed.gddlate.sigma,-2)
    flowertoseed.moisture.mu ~ dnorm(0,0.1)
    flowertoseed.moisture.sigma ~ dunif(0.001,10)
    flowertoseed.moisture.tau <- pow(flowertoseed.moisture.sigma,-2)
    
    
    ## Priors for Year random effects
    for (l in 1:nyears){
    peakflower.year[l] ~dnorm(0,peakflower.tau)
    height.year[l] ~dnorm(0, height.tau)
    }
    
    
    ## Priors for plot random effects
    for (m in 1:nplots){
    peakflower.plot[m] ~dnorm(0,peakflower.tau)
    height.plot[m] ~dnorm(0, height.tau)
    
    }
    
    
    # priors for species-level random effects
    for(k in 1:nspecies){
    height.species[k] ~ dnorm(0, height.tau)
    width.species[k] ~dnorm(0, width.tau)
    peakflower.species[k] ~ dnorm(0, peakflower.tau)
    flowertoseed.species[k] ~dnorm(0, flowertoseed.tau)
    peakflower.SDDDOY.species[k] ~ dnorm(peakflower.SDDDOY.mu, peakflower.SDDDOY.tau)T(0,2)
    peakflower.gddearly.species[k] ~ dnorm(peakflower.gddearly.mu, peakflower.gddearly.tau)
    peakflower.gddlate.species[k] ~ dnorm(peakflower.gddlate.mu,peakflower.gddlate.tau)
    peakflower.moisture.species[k] ~ dnorm(peakflower.moisture.mu,peakflower.moisture.tau)
    flowertoseed.SDDDOY.species[k] ~ dnorm(flowertoseed.SDDDOY.mu, flowertoseed.SDDDOY.tau)
    flowertoseed.gddearly.species[k] ~ dnorm(flowertoseed.gddearly.mu, flowertoseed.gddearly.tau)
    flowertoseed.gddlate.species[k] ~ dnorm(flowertoseed.gddlate.mu,flowertoseed.gddlate.tau)
    flowertoseed.moisture.species[k] ~ dnorm(flowertoseed.moisture.mu,flowertoseed.moisture.tau)
    
    }
    
    # priors for phenological phase parameters
    for(i in 1:nphases){
    width.phase[i] ~ dnorm(0,width.tau)T(-2,5)
    opt.phase[i] ~ dnorm(0,0.1)T(-2.5,2.5)
    }
    
    # likelihood
    for (j in 1:n){
    height[j] <- height.year[year[j]] + height.species[species[j]]
    width[j] <- -1*exp(width.phase[phase[j]] + width.species[species[j]])
    flowertoseed[j] <- flowertoseed.species[species[j]] + (flowertoseed.SDDDOY.species[species[j]] * SDDDOY[j]) + (flowertoseed.gddearly.species[species[j]] *gddearly[j]) +(flowertoseed.gddlate.species[species[j]] *gddlate[j]) + (flowertoseed.moisture.species[species[j]] *moisture[j])
    opt[j] <- peakflower.SDDDOY.species[species[j]] * SDDDOY[j] + (peakflower.gddearly.species[species[j]] *gddearly[j]) +(peakflower.gddlate.species[species[j]] *gddlate[j]) + (peakflower.moisture.species[species[j]] *moisture[j]) + peakflower.plot[plot[j]] + peakflower.year[year[j]] + peakflower.species[species[j]] + 
flowertoseed[j]*step(phase[j] - 1.5)

    logit(y_prob[j]) <- width[j]*(doy[j] - opt[j])^2 + height[j]
    y[j] ~ dbern(y_prob[j])
    }
    
    }
    ",
    file="phenomodel_jags")

#####

### Runs model 

jags_data <- train_data

#Creates a list to move training data from R to JAGS in the right format
jags.data_lags <- list(doy=jags_data$scaledDOY, SDDDOY=jags_data$scaledSDDDOY, gddearly=jags_data$scaledgdd25, gddlate=jags_data$scaledgdd50, moisture=jags_data$scaledsoilmoist, y=jags_data$presence, year=jags_data$year, nyears=length(unique(jags_data$year)), plot=jags_data$plot, nplots=length(unique(jags_data$plot)), phase=jags_data$phase, n=nrow(jags_data),nphases=length(unique(jags_data$phase)), nspecies=length(unique(jags_data$species)), species=jags_data$species)


#Creates list of parameters to monitor
params <- c("peakflower.species", "peakflower.SDDDOY.species","peakflower.gddearly.species", "peakflower.gddlate.species","peakflower.moisture.species","flowertoseed.species", "flowertoseed.SDDDOY.species", "flowertoseed.gddearly.species","flowertoseed.gddlate.species","flowertoseed.moisture.species", "height.year", "height.species", "width.phase", "width.species", "peakflower.plot", "peakflower.year")

#Sets up parallel chains
library(parallel)
cl <- makeCluster(3)


chain3 <- list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 6 )

chain2 <-list(.RNG.name = "base::Mersenne-Twister", .RNG.seed =  4)

chain1 <- list(.RNG.name = "base::Mersenne-Twister", .RNG.seed = 2 )


#Runs model for all data, both flower-fruit and fruit-seed
JAGSresults <- run.jags(model="phenomodel_jags", monitor=params, burnin = 30000, sample=1000, thin=10, inits= list(chain1, chain2, chain3), data=jags.data_lags, n.chains=3, method="bgparallel", cl=cl, keep.jags.files = "JAGSresults")

### Model can be extended if necessary for converging using extend.jags()
