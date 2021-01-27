## WTA Glacier Basin and Skyline Loop testing of ICC using IRR package
## Author: Vicente Velasco
## Date: January 2021
## Description: Includes a function of how to get trailwise
##      icc's for the skyline and glacier basin data

##########################
## 1. Load ICC R package##
##########################

# http://www.cookbook-r.com/Statistical_analysis/Inter-rater_reliability/
install.packages("irr") 	
library(irr)

#######################################
## 2. set directory and read in data ##
#######################################

setwd("C:/Users/vicen/downloads")
rawdata <- read.csv("/Users/vicentevelasco/Downloads/Skyline Loop - Skyline Loop.csv", header = T)
rawComplete <- rawdata[rawdata$pre.2013 == "no",]
rawreviewed <- rawdata[!is.na(rawdata$reviewer.uw),]

n = max(rawdata$orig.position)
m = nrow(rawdata)
orig.position <- rawdata$orig.position

#######################################
## 3. The identifylowirrWTA function ##
#######################################
## This function takes in the WTA trip report data set and the 
## predetermined alpha level of irr that you want to flag 
## and returns a list of the following:
## 1. list of all trail reports with irr's lower than alpha (lowIrrList)
## 2. a list of the all the irr's (allIrrList)
## 3. the mean irr (meanIrr)
## This function assumes that the columns of interest are 20:28

identifylowirrWTA <- function(rawdata, alpha = .8) {
  rawComplete <- rawdata[rawdata$pre.2013 == "no",]
  rawreviewed <- rawdata[!is.na(rawdata$reviewer.uw),]
  
  n = max(rawdata$orig.position)
  m = nrow(rawdata)
  orig.position <- rawdata$orig.position
  lowirrlist <- c()
  allirr <- c()
  for (i in 1:n) {
    #this if statement makes sure that there are two reviewers for the report
    if (sum(rawComplete$orig.position == i) == 2 | sum(rawComplete$orig.position == i) == 3) {
      
      #This if statement tests if the values in the two trail reports are not all NA's 
      if (!all(is.na(rawComplete[rawComplete$orig.position == i, 20:28]))) {
        test <- icc(do.call(rbind, lapply(rawComplete[rawComplete$orig.position == i, 20:28], as.numeric)),model="oneway", type="consistency")
        if (!is.na(test$value)) {
          allirr <- c(allirr, test$value)
          if (test$value < alpha) {
            lowirrlist <- c(lowirrlist, i) 
          }
        }
      }
    }
  }
  meanirr <- mean(allirr)
  lst <- list(lowIrrList = lowirrlist, allIrrList = allirr, meanIrr = meanirr)
  return(lst)
}


lowirrs <- identifylowirrWTA(rawdata, .4)
