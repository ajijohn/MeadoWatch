########################################
## Analysis to calculate flowering season from MeadoWatch Data
## Goals: To quantify the wildflower season from MeadoWatch data, to link with WTA
## This script: 
##     A. Fits curve to each year / plot / species flowering data
##     B. Creates figures 
##     C. (in progress) fit models across species; years
##     D. (in progress) estimate flowering season
#######################################
options(stringsAsFactors = FALSE)

###############################################################
#########Read in Data, assemble important explanatory variables
###############################################################

# Read in data: Phenodat is phenology data; stationdat is lat / long and includes snow disappearance info
PhenoDat <- read.csv("data/MW_PhenoDat_2013_2019.csv", header=TRUE) #phenology data
StationDat <- read.csv("data/MW_SiteDat_2013_2019.csv", header=TRUE) #information about station

# Only the rows in both the files ( Year, Site_Code)
MergePhenoStation <- merge(StationDat,PhenoDat, by=c("Year","Site_Code"))

#Tidy Data into columns needed; order by year
#Tidy Data into columns needed; order by year

# Columns 
colnames(MergePhenoStation)
# [1] "Year"             "Site_Code"        "Transect.x"      
# [4] "Site_Num.x"       "latitude"         "longitude"       
# [7] "elevation"        "meltdate"         "SDD"             
# [10] "Source"           "Transect.y"       "Date"            
# [13] "Month"            "Day"              "Observer"        
# [16] "QA.QC"            "Site_Num.y"       "Species"         
# [19] "Snow"             "Bud"              "Bud_rank"        
# [22] "Flower"           "Flower_rank"      "Fruit"           
# [25] "Fruit_rank"       "Disperse"         "Disperse_rank"   
# [28] "Herb"             "species.notes"    "site.notes"      
# [31] "Data.Entry.Notes"

# Selected  
# "Year"             "Site_Code"        "Transect.x"      
#  "Site_Num.x"       "latitude"         "longitude"
# "SDD" 
# "Date"  
# "Observer"        
#[ "QA.QC"
# "Species"         
# "Snow"             "Bud"     
# "Flower"
# "Fruit" 
# "Disperse"   
# "Herb"
# Selecting relevant columns - see above
PhenoSite_0 <- MergePhenoStation[,c(1:6, 9, 12, 15:16, 18:20, 22, 24, 26, 28)]
# New columns are
# [1] "Year"       "Site_Code"  "Transect.x" "Site_Num.x" "latitude"  
# [6] "longitude"  "SDD"        "Date"       "Observer"   "QA.QC"     
# [11] "Species"    "Snow"       "Bud"        "Flower"     "Fruit"     
# [16] "Disperse"   "Herb"   
# Reorder by Year
PhenoSite_0 <- PhenoSite_0[order(PhenoSite_0[,1]),]

#Calculate Julian Days of observations; DSS=days since snow; add to PhenoSite
# Picking up unique years
# 7 years - 2013 - 2019
Yrs <- unique(PhenoSite_0$Year);
# DOY - Day of the Year (DOY) assign to empty vector
DOY <- c()

# Convert the 'Date' - observed date to Julian (DOY)
for(i in 1:length(Yrs)){
  tmpdates <- PhenoSite_0$Date[PhenoSite_0$Year==Yrs[i]]
  Jan1Jul <- as.Date(paste(Yrs[i],"-01-01", sep=""))
  ObsJulDayYr <- julian(as.Date(as.character(tmpdates),"%m/%d/%Y"), origin=Jan1Jul)
  DOY <- c(DOY,ObsJulDayYr)
}

# Adding the DOY column to merged reordered dataframe 
PhenoSite_0 <- cbind(PhenoSite_0, DOY)
# reorder from 18 columns
# [1] "Year"       "Site_Code"  "Transect.x" "Site_Num.x" "latitude"  
# [6] "longitude"  "SDD"        "Date"       "Observer"   "QA.QC"     
# [11] "Species"    "Snow"       "Bud"        "Flower"     "Fruit"     
# [16] "Disperse"   "Herb"       "DOY"  

# to 10 columns 
#[1] "Year"       "Transect.x" "Site_Code"  "Site_Num.x" "Date"      
#[6] "DOY"        "SDD"        "QA.QC"      "Species"    "Flower" 
PhenoSite <- PhenoSite_0[,c(1,3,2,4,8,18,7,10,11,14)] #reorganize data


###############################################################
#########Define Maximum Likelihood Models (to fit to data)
###############################################################

##INDIVIDUAL SPECIES MODELS
#Null model - assumes the probability of flowering does not vary with time
nullfit <- function (param){ 
	meanp  <- param[1]
	pred   <- rep(meanp, times=length(phenophase))
	llik     <- dbinom(phenophase,1,pred, log=TRUE)
	return(-sum(llik))
}

#Alternative model - assumes flowering varies with DOY according to a guassian curve
#can fit this to individual plot / year / species data, data for a species in all plots, etc
curvefit <- function (param){ #curve fitting function for mle
	 peakp  <- param[1]
	 rangep <- param[2]
	 maxp   <- param[3]
	 pred   <- maxp*(1/(rangep*sqrt(2*pi)))*exp(-0.5*((days-peakp)/rangep)^2)
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
  pred   <- maxp*(1/(rangep*sqrt(2*pi)))*exp(-0.5*((days-peakp)/rangep)^2)
  llik     <- dbinom(phenophase,1,pred, log=TRUE)
  return(-sum(llik))
}

#Fits data to all species within a trail & a season, plot & species effects
#haven't tried this yet - needs trouble shooting
curvefit_plotspmod <- function (param){ #curve fitting function for mle
  peakp_ref  <- param[1:nspp] #species specific peak in reference plot
  rangep_all <- param[(nspp+1):(nspp*2)] #species specific range in all plots
  maxp_all   <- param[(2*nspp+1):(nspp*3)] #species maxp in all plots
  peakp_plot <- param[(3*nspp+1):(3*nspp+nplt)] # plot effects, assume no interaction
  maxp <- maxp_all[spp]
  rangep <- rangep_all[spp]
  peakp <- peakp_ref[spp]
  peakp_add <- peakp_plot[plt]; peakp_add[is.na(peakp_add)] <- 0 #add plot level effects
  peakp <- peakp + peakp_add
  
  pred   <- maxp*(1/(rangep*sqrt(2*pi)))*exp(-0.5*((days-peakp)/rangep)^2)
  llik     <- dbinom(phenophase,1,pred, log=TRUE)
  return(-sum(llik))
}



################################################################################
#########Define function to plot phenological curve given optim, range, maxp ####
#################################################################################


#Predicts probability of observing phenological phase, given parameters (used for plotting)
predphen <- function (xx, param){ 
  days    <- xx
	peakp  <- param[1]
	rangep <- param[2]
	maxp   <- param[3]
	pred   <- maxp*(1/(rangep*sqrt(2*pi)))*exp(-0.5*((days-peakp)/rangep)^2)
	return(pred)
}


#############################################################
######### Fit model per plot / species year  ##############
############################################################

#Define species of interest (for reference, RL and GB focal species listed)
species <- c("ANOC","ARLA","ASLE","CAMI","CAPA","ERGR","ERMO","ERPE","LUAR","PEBR","VASI"); nspp <- length(species)
#REFLECTION LAKES SPECIES ARE <- c("ANOC","CAPA","ERMO","ERPE","LIGR","LUAR","MIAL","PEBR","POBI","VASI")
#GLACIER BASIN SPECIES ARE <- c("ANAR","ARLA","ASLE","CAMI", "ERGR","LUAR","MEPA","PEBR","POBI","VASI")
#Exclude non focal species of interest
specieskeep <- PhenoSite$Species %in% species
PhenoSite_Focal <- PhenoSite[specieskeep,]

#how many years?
years <- unique(PhenoSite_Focal$Year)

#where to save output & SDD of trail-species-plot-year specific fits
parameters <- c() 

##Now nested for loops to fit curves for all years, focal species, plots within trails
for(i in 1:length(years)){ #loop for each year
  #extract data for that year
  PhenoSite_Year <- PhenoSite_Focal[PhenoSite_Focal$Year==years[i],]
  
  #pull out unique plots for that year; varies by year
  plots <- unique(PhenoSite_Year$Site_Code)
  
  #now fit data per plot, per species
  for(j in 1:length(plots)){ # loop for each plot
    PhenoSite_YearPlot <- PhenoSite_Year[PhenoSite_Year$Site_Code==plots[j],]

    #Identify species in plot
    speciesinplot <- unique(PhenoSite_YearPlot$Species)
    if(length(speciesinplot)==0){next} #break if focal spp not in plot
    
    #loop for each species in the plot
    for(k in 1:length(speciesinplot)){

      #pull out data for species in plot      
      PhenoSite_YearPlotSpecies <- PhenoSite_YearPlot[PhenoSite_YearPlot$Species==speciesinplot[k],]
      
      #define parameters for curvefitting
      days <- PhenoSite_YearPlotSpecies$DOY #explanatory variable: DOY 
      phenophase <- PhenoSite_YearPlotSpecies$Flower #Response variables: flowers

      #If less than 5 total observations, do not fit
      if(length(days)<6){next}
      #remove days when no observations were made; NA in phenophase
	    days <- days[is.na(phenophase)==FALSE]
	    phenophase <- phenophase[is.na(phenophase)==FALSE]
	    
	    #add three weeks of zeroes before earliest SDD in those plots
	    SDDplt <- min(PhenoSite_YearPlotSpecies$SDD)
	    days <- c(SDDplt-21, SDDplt-14, SDDplt-7,days)
	    phenophase <- c(0,0,0,phenophase)

      #now fit null
      model0 <- optimize(nullfit, c(0.000001,0.999999)) #fit null model
	
      #now fit alternative model - curve
      param <- c(mean(days), sd(days), 0.5) # initial parameters: model fits pretty fussy about this
      model1 <- optim(param, curvefit, control = list(maxit = 50000))
      if(model1$convergence==1){print(paste(speciesinplot[k],"no convergence", sep="-"))}
    	
      #calculate AIC, p value
      AICnull <- round(2*(model0$objective+1),1)
      AICalt <- round(2*(model1$value + 3),1)
      pcurve <- signif(pchisq(model1$value-model0$objective,2),3)
    
      #print output if desired
      #print(speciesinplot[k])
      #print(paste("AICnull=",AICnull,"AICalt=",AICalt,"p(curve)=",pcurve))
	
      #save parameters to plot all species in plot
      tmp_pars <- c(years[i], as.character(plots[j]), SDDplt, as.character(speciesinplot[k]), model1$par[1:3])
      parameters <- rbind(parameters, tmp_pars)
    }
  }
}
  
#turn parameters into a data frame

# Range - Width
# Max - Height
# Peak - Day of the year of Peak flowerimg

dimnames(parameters) <- list(c(), c("year","plot","SDD","species","peak","range","max"))
parameters <- data.frame(parameters)

#change storage type to numeric - all except plot (since a few plots have a, b)
parameters$year <- as.numeric(parameters$year)
parameters$species <- as.factor(parameters$species)
parameters$SDD <- as.numeric(parameters$SDD)
parameters$peak <- as.numeric(parameters$peak)
parameters$range <- as.numeric(parameters$range)
parameters$max <- as.numeric(parameters$max)


#############################################################################
####  Create plots from plot/species fits; runs some simple stats  ##########
#############################################################################

###############
##Plot all curves for all species and all plots per year and per trail

#define plotting colors, sufficient for all focal species
plotcol <- c("pink","orangered","yellow","purple","lightblue","grey","magenta","yellowgreen",
             "navyblue","azure4","yellow4","yellowgreen","orchid","turquoise","salmon","maroon")
for(trail in 1:2){
  if(trail==1){parameters2 <- parameters[substr(parameters$plot,1,2)=="RL",]}
  if(trail==2){parameters2 <- parameters[substr(parameters$plot,1,2)=="GB",]}
  
  #how many years? differs per trail
  years <- unique(parameters2$year)
  
  X11(width=4.5, height=9)
  par(mfrow=c(length(years),1),omi=c(0,0,0,0), mai=c(0.1,0.3,0.1,0.0),tck=-0.01, mgp=c(1.25,0.25,0),xpd=TRUE)

  #plot per year
  for(i in 1:length(years)){
    paryear <- parameters2[parameters2$year==years[i],]
  
    plot(185,0.5, xlim=c(105,285), ylim=c(0,1.01),type="n",
         xaxp=c(150,270,5), xaxt="n", xlab="Time",ylab="Flowering")
    text(c(105,135,165,195,225,255),-0.1, labels=c("Apr","May","Jun","Jul","Aug","Sep"))
    text(110,0.9,labels=years[i])
    if(i==1&trail==1){title("Reflection Lakes")}
    if(i==1&trail==2){title("Glacier Basin")}
    
    #now pull out all data for a species
    for(j in 1:length(species)){
      paryearsp <- paryear[paryear$species==species[j],]
      if(dim(paryearsp)[1]==0){next}

      #now pull out data per species, and plot
      for(k in 1:dim(paryearsp)[1]){
         xx <- seq(105,285)
         n <- which(dimnames(paryearsp)[[2]]=="peak")
         pars <- unlist(paryearsp[k,n:(n+2)])
         yy <- predphen(xx, pars)
         if(max(yy)>1){yy <- yy/(max(yy))} #bit of a kluge - force max to be < 1
         xx2 <- xx[yy[]>=0.001]; yy2 <- yy[yy[]>=0.001]
         lines(xx2,yy2, col=plotcol[j],lwd=2)
      }
    }
  }
}

#######################
##Plot peak flowering estimates of all species by year, both trails on 1 plot
X11(width=7.5, height=4)
par(mfrow=c(1,1),omi=c(0,0,0,0), mai=c(0.5,0.4,0.4,0.2),tck=-0.01, mgp=c(1.25,0.25,0),xpd=TRUE)

#create plot to add points to
earlypk <- min(parameters$peak); latepk <- max(parameters$peak)
plot(2016,175, xlim=c(2012, 2020), ylim=c(earlypk,latepk),type="n",
 xaxp=c(2012,2020,8), yaxt="n", xlab="Year",ylab="Flowering")
text(2011.5,srt=90, c(135,165,195,225),-0.1, labels=c("May","Jun","Jul","Aug"))
  

for(trail in 1:2){
  if(trail==1){parameters2 <- parameters[substr(parameters$plot,1,2)=="RL",]}
  if(trail==2){parameters2 <- parameters[substr(parameters$plot,1,2)=="GB",]}
  
  #how many years? differs per trail
  years <- unique(parameters2$year)

  #extract data per year, plot
  for(i in 1:length(years)){
    paryear <- parameters2[parameters2$year==years[i],]
 
    #now pull out all data for a species
    for(j in 1:length(species)){
      paryearsp <- paryear[paryear$species==species[j],]
      if(dim(paryearsp)[1]==0){next}
      pks <- paryearsp$peak
      tiny <- 0; if(trail==2){tiny <- 0.333}
      if(trail==1){pltshp <-21}; if(trail==2){pltshp <- 24}
      points((rep(years[i],length(pks))+jitter(rep(tiny,length(pks)),2.5)),pks,pch=pltshp, bg=plotcol[j], cex=1.25)
    }
  }
}

legend(x="topleft", legend=c("RL","GB"), pch=c(21,24), pt.bg="gray",cex=0.75)


##One panel per year, one graph per plot per trail; all species plotted on that graph
#plotting parameters
years <- unique(parameters$year)

#one plot per year, a
for(i in 1:length(years)){
  paryear <- parameters[parameters$year==years[i],]
  plts <-unique(paryear$plot); nplts <- length(plts)
  
  X11(width=9, height=7)
  if(nplts<10){par(mfrow=c(3,3),omi=c(0,0,0,0), mai=c(0.1,0.1,0.1,0.0),tck=-0.01, mgp=c(1.25,0.25,0),xpd=TRUE)}
  if(nplts>9 & nplts<26){par(mfrow=c(5,5),omi=c(0,0,0,0), mai=c(0.1,0.1,0.1,0.0),tck=-0.01, mgp=c(1.25,0.25,0),xpd=TRUE)}
  if(nplts>25){par(mfrow=c(6,5),omi=c(0,0,0,0), mai=c(0.1,0.1,0.1,0.0),tck=-0.01, mgp=c(1.25,0.25,0),xpd=TRUE)}
  
  #now pull out all data for a plot
  for(j in 1:nplts){
    paryearplt <- paryear[paryear$plot==plts[j],]
    splot <- unique(paryearplt$species); nsp <- length(splot)
    
    plot(185,0.5, xlim=c(135,285), ylim=c(0,1.01),type="n",
       xaxp=c(150,240,4), xaxt="n", xlab="Time",ylab="Flowering")
    text(c(135,165,195,225,255),-0.075, labels=c("May","Jun","Jul","Aug","Sept"))
    text(135,0.95,adj=0, labels=paste(years[i],"plot",plts[j],sep=" "))
  
    #now plot
    for(k in 1:nsp){
      xx <- seq(135,285)
      n <- which(dimnames(paryearplt)[[2]]=="peak")
      pars <- unlist(paryearplt[k,n:(n+2)])
      yy <- predphen(xx, pars)
      if(max(yy)>1){yy <- yy/(max(yy))} #bit of a kluge - force max to be < 1
      xx2 <- xx[yy[]>=0.001]; yy2 <- yy[yy[]>=0.001]
      whichsp <- which(species[]==splot[k])
      lines(xx2,yy2, col=plotcol[whichsp],lwd=2)
    }
  }
}


##One graphics window per year; all species specific plot curves on one graph
#plotting parameters

#plot per year
for(i in 1:length(years)){
  paryear <- parameters[parameters$year==years[i],]

  X11(width=7, height=7)
  if(length(species)<=9){par(mfrow=c(3,3),omi=c(0,0,0,0), mai=c(0.2,0.3,0.3,0.0),tck=-0.01, mgp=c(1.25,0.25,0),xpd=TRUE)}
  if(length(species)>9){par(mfrow=c(4,3),omi=c(0,0,0,0), mai=c(0.2,0.3,0.3,0.0),tck=-0.01, mgp=c(1.25,0.25,0),xpd=TRUE)}
  
          
  #now pull out all data for a species
  for(j in 1:length(species)){
    paryearsp <- paryear[paryear$species==species[j],]
    if(dim(paryearsp)[1]==0){next}
    
    plot(185,0.5, xlim=c(135,285), ylim=c(0,1.01),type="n",
       xaxp=c(150,240,4), xaxt="n", xlab="Time",ylab="Flowering")
    text(c(135,165,195,225,255),-0.075, labels=c("May","Jun","Jul","Aug","Sep"))
    title(paste(years[i],species[j], sep=" "))
 
    #now plot
    for(k in 1:dim(paryearsp)[1]){
      xx <- seq(135,285)
      n <- which(dimnames(paryearsp)[[2]]=="peak")
      pars <- unlist(paryearsp[k,n:(n+2)])
      yy <- predphen(xx, pars)
      if(max(yy)>1){yy <- yy/(max(yy))} #bit of a kluge - force max to be < 1
      xx2 <- xx[yy[]>=0.001]; yy2 <- yy[yy[]>=0.001]
      lines(xx2,yy2, col=plotcol[j],lwd=2)
    }
  }
}

##Assess how parameters vary with SDD, trail and year
#In same for loop, create graphs (one per species; parameter); run some tests

yearcol <- c("purple","grey","orange","lightblue","yellowgreen","pink","navyblue") #colors per year
testparsAIC <- c() #output of simple exploratory models

for(i in 1:length(species)){
  parspecies <- parameters[parameters$species==species[i],]
  
  #set shape
  pltshp <- rep(21, length=dim(parspecies)[1])
  pltshp[substr(parspecies$plot,1,2)=="GB"] <- 24
  
  #setup plot
  X11(width=8.5, height=4.0)
  par(mfrow=c(1,3),omi=c(0,0,0,0), mai=c(0.35,0.335,0.3,0.0),tck=-0.02, mgp=c(1.25,0.25,0),xpd=TRUE)
  
  #SDD vs peak
  plot(parspecies$SDD,parspecies$peak, xlab="SDD",ylab="peak", pch=pltshp, bg=yearcol[parspecies$year-2012], cex=1.5)
  title(paste(species[i],"- SDD vs. optim"))
  legend(x="bottomright",c("'13","'14","'15","'16","'17","'18","'19","GB","RL"), 
         pch=c(21,21,21,21,21,21,21,24,21), pt.bg=c(yearcol,"gray","gray"))

  #test for plotting - how strongly are peak and SDD correlated?
  optimtest <- cor.test(parspecies$SDD,parspecies$peak)
  mtext(paste("cor=",round(optimtest$estimate,3),sep=""), side=3, adj=0.025, cex=0.75,line=-1)
  mtext(paste("p=",round(optimtest$p.value,3),sep=""), side=3, adj=0.025, cex=0.75,line=-2)
  
  #test for output - does peak flowering depend on SDD, year, trail?
  ntrails <- length(unique(substr(parspecies$plot,1,2)))
  if(ntrails==1){
    testall <- lm(parspecies$peak ~ parspecies$SDD + parspecies$year)
    testSDD <- lm(parspecies$peak ~ parspecies$SDD)
    testyear <- lm(parspecies$peak ~ parspecies$year)
    testnull <- lm(parspecies$peak ~ 1)
    tmpAIC <- c(AIC(testnull), AIC(testyear),AIC(testSDD),NA, AIC(testall))
    }
  if(ntrails==2){
    testall <- lm(parspecies$peak ~ parspecies$SDD + parspecies$year + as.factor(substr(parspecies$plot,1,2)))
    testtrail <- lm(parspecies$peak ~ as.factor(substr(parspecies$plot,1,2)))
    testSDD <- lm(parspecies$peak ~ parspecies$SDD)
    testyear <- lm(parspecies$peak ~ parspecies$year)
    testnull <- lm(parspecies$peak ~ 1)
    tmpAIC <- c(AIC(testnull), AIC(testyear),AIC(testSDD),AIC(testtrail), AIC(testall))
  }

  #save AIC
  tmpoutput <- c(species[i],"peak", tmpAIC)
  testparsAIC <- rbind(testparsAIC,tmpoutput)
    
  #SDD vs duration / range  
  plot(parspecies$SDD,parspecies$range, xlab="SDD",ylab="range", pch=pltshp, bg=yearcol[parspecies$year-2012], cex=1.5)
  title(paste(species[i],"- SDD vs. range"))
  
  #test - how strongly are range and SDD correlated?
  optimtest <- cor.test(parspecies$SDD,parspecies$range)
  mtext(paste("cor=",round(optimtest$estimate,3),sep=""), side=3, adj=0.025, cex=0.75,line=-1)
  mtext(paste("p=",round(optimtest$p.value,3),sep=""), side=3, adj=0.025, cex=0.75,line=-2)
  
  #test for output - does duration flowering depend on SDD, year, trail?
  if(ntrails==1){
    testall <- lm(parspecies$range ~ parspecies$SDD + parspecies$year)
    testSDD <- lm(parspecies$range ~ parspecies$SDD)
    testyear <- lm(parspecies$range ~ parspecies$year)
    testnull <- lm(parspecies$range ~ 1)
    tmpAIC <- c(AIC(testnull), AIC(testyear),AIC(testSDD),NA, AIC(testall))
  }
  if(ntrails==2){
    testall <- lm(parspecies$range ~ parspecies$SDD + parspecies$year + as.factor(substr(parspecies$plot,1,2)))
    testtrail <- lm(parspecies$range ~ as.factor(substr(parspecies$plot,1,2)))
    testSDD <- lm(parspecies$range ~ parspecies$SDD)
    testyear <- lm(parspecies$range ~ parspecies$year)
    testnull <- lm(parspecies$range ~ 1)
    tmpAIC <- c(AIC(testnull), AIC(testyear),AIC(testSDD),AIC(testtrail), AIC(testall))
  }
  
  #save AIC
  tmpoutput <- c(species[i],"duration", tmpAIC)
  testparsAIC <- rbind(testparsAIC,tmpoutput)

  #SDD vs max height parameter
  plot(parspecies$SDD,parspecies$max, xlab="SDD",ylab="max", pch=pltshp, bg=yearcol[parspecies$year-2012], cex=1.5)
  title(paste(species[i],"- SDD vs. max height"))
  
  #test - how strongly are range and SDD correlated?
  optimtest <- cor.test(parspecies$SDD,parspecies$max)
  mtext(paste("cor=",round(optimtest$estimate,3),sep=""), side=3, adj=0.025, cex=0.75,line=-1)
  mtext(paste("p=",round(optimtest$p.value,3),sep=""), side=3, adj=0.025, cex=0.75,line=-2)

  #test for output - does max flowering depend on SDD, year, trail?
  if(ntrails==1){
    testall <- lm(parspecies$max ~ parspecies$SDD + parspecies$year)
    testSDD <- lm(parspecies$max ~ parspecies$SDD)
    testyear <- lm(parspecies$max ~ parspecies$year)
    testnull <- lm(parspecies$max ~ 1)
    tmpAIC <- c(AIC(testnull), AIC(testyear),AIC(testSDD),NA, AIC(testall))
  }
  if(ntrails==2){
    testall <- lm(parspecies$max ~ parspecies$SDD + parspecies$year + as.factor(substr(parspecies$plot,1,2)))
    testtrail <- lm(parspecies$max ~ as.factor(substr(parspecies$plot,1,2)))
    testSDD <- lm(parspecies$max ~ parspecies$SDD)
    testyear <- lm(parspecies$max ~ parspecies$year)
    testnull <- lm(parspecies$max ~ 1)
    tmpAIC <- c(AIC(testnull), AIC(testyear),AIC(testSDD),AIC(testtrail), AIC(testall))
  }
  
  #save AIC
  tmpoutput <- c(species[i],"max", tmpAIC)
  testparsAIC <- rbind(testparsAIC,tmpoutput)
}

#maketestparsAIC a data frame
dimnames(testparsAIC) <- list(c(), c("species","parameter","AICnull","AICyear","AICSDD","AICtrail","AICall"))
testparsAIC <- data.frame(testparsAIC)

#change storage type to numeric - all except plot (since a few plots have a, b)
testparsAIC$species <- as.factor(testparsAIC$species)
testparsAIC$parameter <- as.factor(testparsAIC$parameter)
testparsAIC$AICnull <- as.numeric(testparsAIC$AICnull)
testparsAIC$AICyear <- as.numeric(testparsAIC$AICyear)
testparsAIC$AICtrail <- as.numeric(testparsAIC$AICtrail)
testparsAIC$AICSDD <- as.numeric(testparsAIC$AICSDD)
testparsAIC$AICall <- as.numeric(testparsAIC$AICall)


  

##########Now fit model where peak flowering is based on SDD
species <- c("ANOC","ARLA","ASLE","CAMI","CAPA","ERGR","ERMO","ERPE","LUAR","PEBR","VASI"); nspp <- length(species)
#REFLECTION LAKES SPECIES ARE <- c("ANOC","CAPA","ERMO","ERPE","LIGR","LUAR","MIAL","PEBR","POBI","VASI")
#GLACIER BASIN SPECIES ARE <- c("ANAR","ARLA","ASLE","CAMI", "ERGR","LUAR","MEPA","PEBR","POBI","VASI")
#Exclude non focal species of interest
specieskeep <- PhenoSite$Species %in% species
PhenoSite_Focal <- PhenoSite[specieskeep,]

#where to save output
parameters_snowmod <- c() #save parameters for per species fit, peak f based on SDD

##Now nested for loops to fit curves for all years, species, plots
for(i in 1:length(species)){ #loop for each year
  #extract data for that year
  PhenoSite_Species <- PhenoSite_Focal[PhenoSite_Focal$Species==species[i],]
  
  #assemble data - need to add zeros before SDD, define trail year combos
  yrplt <- paste(PhenoSite_Species$Year, PhenoSite_Species$Site_Code, sep="")
  yrplts <- unique(yrplt)
  yeartrail <- unique(substr(yrplts,1,6)); yeartrail <- yeartrail[order(yeartrail)]
  ntrlyr <- length(yeartrail)
  trlyr <- c()
  days <- c()
  phenophase <- c()
  SDD <- c()
  
  for(j in 1:length(yrplts)){
     PhenoSite_YearPlotSpecies <- PhenoSite_Species[yrplt[]==yrplts[j],]
     tmpdays <- PhenoSite_YearPlotSpecies$DOY #explanatory variable: DOY 
     tmpphenophase <- PhenoSite_YearPlotSpecies$Flower #Response variables: flowers
     tmpSDD <- PhenoSite_YearPlotSpecies$SDD
     if(length(tmpdays)<10){next}
     #remove days when no observations were made; NA in phenophase
     tmpdays <- tmpdays[is.na(tmpphenophase)==FALSE]
     tmpSDD <- tmpSDD[is.na(tmpphenophase)==FALSE]
     tmpphenophase <- tmpphenophase[is.na(tmpphenophase)==FALSE]
     SDDplt <- min(tmpSDD)
     tmpdays <- c(SDDplt-21, SDDplt-14, SDDplt-7,tmpdays)
     tmpphenophase <- c(0,0,0,tmpphenophase)
     tmpSDD <- c(tmpSDD[1:3],tmpSDD)
    
     #determine which trail year combo
     whichtrlyr <- substr(paste(PhenoSite_YearPlotSpecies$Year[1], PhenoSite_YearPlotSpecies$Site_Code[1], sep=""),1,6)
     tmptrlyr <- rep(which(whichtrlyr==yeartrail[]),times=length(tmpdays))
     #now append to days, phenophase,trlyear
     days <- c(days, tmpdays)
     phenophase <- c(phenophase, tmpphenophase)
     SDD <- c(SDD, tmpSDD)
     trlyr <- c(trlyr, tmptrlyr)
  }
  #now fit alternative model - optim as function of SDD
  #first use parameters to come up with reasonable priors
  #slope & intercept for SDD vs optim
  
  optimcoef <- unlist(coef(lm(parameters$peak~parameters$SDD)))
  sprng <- mean(parameters$range[parameters$species==species[i]])
  spmx <- mean(parameters$max[parameters$species==species[i]])

  #set parameters & fit model
  param <- c(optimcoef, rep(sprng, times=ntrlyr), rep(spmx, times=ntrlyr)) # initial parameters: model fits pretty fussy about this
  modelsnowmod <- optim(param, curvefit_snowmod, control = list(maxit = 50000))

  #Determine which years / trails (for plotting, saving)
  trails <- "both"
  if(length(unique(substr(yeartrail,5,6)))==1){
    if(unique(substr(yeartrail,5,6))=="RL"){trails <- "Reflection Lakes"}
    if(unique(substr(yeartrail,5,6))=="GB"){trails <- "Glacier Basin"}
  }
  #print(species[i]); print(modelsnowmod); print(modelsnowmod$convergence)
  ifelse(trails=="Glacier Basin", years <- seq(2015,2019), years <- seq(2013,2019))
  
  #Save parameters to 
  
}

#Fits model per species using snow disappearance as prediction - across trails / plots / years
#haven't tried this yet - needs trouble shooting
curvefit_snowmod <- function (param){ #curve fitting function for mle
  intpeak  <- param[1]
  slopepeak <- param[2]
  rangep_all <- param[3:(2+ntrlyr)] #different range per trail, year
  maxp_all   <- param[(3+ntrlyr):(2+2*ntrlyr)] #different max per trail, year
  peakp <- intpeak + slopepeak*SDD
  rangep <- rangep_all[trlyr]
  maxp <- maxp_all[trlyr]
  pred   <- maxp*(1/(rangep*sqrt(2*pi)))*exp(-0.5*((days-peakp)/rangep)^2)
  llik     <- dbinom(phenophase,1,pred, log=TRUE)
  return(-sum(llik))
}


      
