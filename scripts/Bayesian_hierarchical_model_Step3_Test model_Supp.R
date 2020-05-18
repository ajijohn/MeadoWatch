## Author: Meera Lee Sethi
## Last worked on: 12-27-2019
## Script for testing the Bayesian hierarchical model used in "Early snowmelt and warmer, drier summers shrink post-flowering transition times in subalpine wildflowers"

#Sets up workspace

library(PRROC)
setwd("~/Google Drive/Research/Wildflower Reproductive Phenology")

testdata <- read.csv("./final_test_data.csv") #Reads in held-back data for testing
summary <- read.csv("./coefficients.csv") #Reads in file containing model fit coefficients


antilogit <- function(x){exp(x)/(1 + exp(x))} #function to convert data points from a scale of -infinity to +infinity (the scale of the linear predictor) to a scale of 0-1 (the probability scale)

## testing using final median values

preds <- matrix(ncol=3, nrow = 5054)

  for(j in 1:5054){
    year = testdata[j,24]
    species = testdata[j,25]
    phase = testdata[j,22]
    plot = testdata[j,23]
    SDDDOY= testdata[j,27]
    DOY = testdata[j,26]
    gddearly= testdata[j,28]
    gddlate= testdata[j,29]
    moisture= testdata[j,30]
    presence = testdata[j,18]

    height.year = summary$Median[summary$Variables==paste("height.year[", year, "]", sep = "")]
    height.species =summary$Median[summary$Variables == paste("height.species[", species, "]", sep = "")]
    height = height.year + height.species
    width.species = summary$Median[summary$Variables==paste("width.species[", species, "]", sep = "")]
    width.phase = summary$Median[summary$Variables == paste("width.phase[", phase, "]", sep = "")]
    width = -1*exp(width.species + width.phase)
    peakflower.species = summary$Median[summary$Variables == paste("peakflower.species[", species, "]", sep = "")]
    peakflower.year = summary$Median[summary$Variables== paste("peakflower.year[", year, "]", sep = "")]
    peakflower.species = summary$Median[summary$Variables== paste("peakflower.species[", species, "]", sep = "")]
    peakflower.plot = summary$Median[summary$Variables == paste("peakflower.plot[", plot, "]", sep = "")]
    peakflower.SDDDOY.species = summary$Median[summary$Variables == paste("peakflower.SDDDOY.species[", species, "]", sep = "")]
    flowertoseed.species = summary$Median[summary$Variables == paste("flowertoseed.species[", species, "]", sep = "")]
    flowertoseed.SDDDOY.species = summary$Median[summary$Variables == paste("flowertoseed.SDDDOY.species[", species, "]", sep = "")]
    flowertoseed.gddearly.species = summary$Median[summary$Variables== paste("flowertoseed.gddearly.species[", species, "]", sep = "")]
    flowertoseed.gddlate.species = summary$Median[summary$Variables == paste("flowertoseed.gddlate.species[", species, "]", sep = "")]
    flowertoseed.moisture.species =  summary$Median[summary$Variables== paste("flowertoseed.moisture.species[", species, "]", sep = "")]
    flowertoseed = flowertoseed.species+ (flowertoseed.SDDDOY.species * SDDDOY) + (flowertoseed.gddearly.species *gddearly) +(flowertoseed.gddlate.species *gddlate) + (flowertoseed.moisture.species *moisture)
    opt = peakflower.SDDDOY.species*SDDDOY + peakflower.plot + peakflower.year + peakflower.species +
      ifelse(phase>=2,1,0)*flowertoseed  #if phenophase is 2 or higher, add first lag
    y = width*(DOY- opt)^2 + height
    y_prob = antilogit(y)
    y.pred = rbinom(1, size = 1, prob = y_prob)
    preds[j,2] <- presence
    preds[j,3] <-y_prob
    preds[j,1] <- y.pred
  }

hist(preds[,3][preds[,2]==1]) # model probabilities when phase is present
hist(preds[,3][preds[,2]==0]) # model probabilities when phase is absent

#write.csv(file="predictions.csv", x=preds)

#preds <- as.data.frame(read.csv(file="predictions.csv"))


# ROC Curve    
roc <- roc.curve(scores.class0 = fg, scores.class1 = bg, curve = T)
plot(roc)

# PR Curve
pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)

plot(pr)


## testing for overfitting

preds_training <- matrix(ncol=3, nrow = 45430)

for(j in 1:45430){
  year = traindata[j,24]
  species = traindata[j,25]
  phase = traindata[j,22]
  plot = traindata[j,23]
  SDDDOY= traindata[j,27]
  DOY = traindata[j,26]
  gddearly= traindata[j,28]
  gddlate= traindata[j,29]
  moisture= traindata[j,30]
  presence = traindata[j,18]
  
  height.year = summary$Median[summary$Variables==paste("height.year[", year, "]", sep = "")]
  height.species =summary$Median[summary$Variables == paste("height.species[", species, "]", sep = "")]
  height = height.year + height.species
  width.species = summary$Median[summary$Variables==paste("width.species[", species, "]", sep = "")]
  width.phase = summary$Median[summary$Variables == paste("width.phase[", phase, "]", sep = "")]
  width = -1*exp(width.species + width.phase)
  peakflower.species = summary$Median[summary$Variables == paste("peakflower.species[", species, "]", sep = "")]
  peakflower.year = summary$Median[summary$Variables== paste("peakflower.year[", year, "]", sep = "")]
  peakflower.species = summary$Median[summary$Variables== paste("peakflower.species[", species, "]", sep = "")]
  peakflower.plot = summary$Median[summary$Variables == paste("peakflower.plot[", plot, "]", sep = "")]
  peakflower.SDDDOY.species = summary$Median[summary$Variables == paste("peakflower.SDDDOY.species[", species, "]", sep = "")]
  flowertoseed.species = summary$Median[summary$Variables == paste("flowertoseed.species[", species, "]", sep = "")]
  flowertoseed.SDDDOY.species = summary$Median[summary$Variables == paste("flowertoseed.SDDDOY.species[", species, "]", sep = "")]
  flowertoseed.gddearly.species = summary$Median[summary$Variables== paste("flowertoseed.gddearly.species[", species, "]", sep = "")]
  flowertoseed.gddlate.species = summary$Median[summary$Variables == paste("flowertoseed.gddlate.species[", species, "]", sep = "")]
  flowertoseed.moisture.species =  summary$Median[summary$Variables== paste("flowertoseed.moisture.species[", species, "]", sep = "")]
  flowertoseed = flowertoseed.species+ (flowertoseed.SDDDOY.species * SDDDOY) + (flowertoseed.gddearly.species *gddearly) +(flowertoseed.gddlate.species *gddlate) + (flowertoseed.moisture.species *moisture)
  opt = peakflower.SDDDOY.species*SDDDOY + peakflower.plot + peakflower.year + peakflower.species +
    ifelse(phase>=2,1,0)*flowertoseed  #if phenophase is 2 or higher, add first lag
  y = width*(DOY- opt)^2 + height
  y_prob = antilogit(y)
  y.pred = rbinom(1, size = 1, prob = y_prob)
  preds_training[j,2] <- presence
  preds_training[j,3] <-y_prob
  preds_training[j,1] <- y.pred
}

fg <- preds_training[,3][traindata$presence == 1]
bg <- preds_training[,3][traindata$presence == 0]

# PR Curve
plot(pr)
