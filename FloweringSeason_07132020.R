########################################
## Analysis to calculate flowering season from MeadoWatch Data
## Goals: To quantify the wildflower season from MeadoWatch data, to link with WTA
## This script: 
##     A. Fits curve to each year / plot / species flowering data
##     B. Creates figures 
##     C. (in progress) fit models across species; years
##     D. (in progress) estimate flowering season
## Last worked on by Janneke Hille Ris Lambers: July 13, 2020
#######################################


###############################################################
#########Read in Data, assemble important explanatory variables
###############################################################

# Read in data: Phenodat is phenology data; stationdat is lat / long and includes snow disappearance info
PhenoDat <- read.csv("data/MW_PhenoDat_2013_2019.csv", header=TRUE) #phenology data
StationDat <- read.csv("data/MW_SiteDat_2013_2019.csv", header=TRUE) #information about station
MergePhenoStation <- merge(StationDat,PhenoDat, by=c("Year","Site_Code"))

#Tidy Data into columns needed; order by year
PhenoSite_0 <- MergePhenoStation[,c(1:6, 9, 12, 15:16, 18:20, 22, 24, 26, 28)]
PhenoSite_0 <- PhenoSite_0[order(PhenoSite_0[,1]),]

#Calculate Julian Days of observations; DSS=days since snow; add to PhenoSite
Yrs <- unique(PhenoSite_0$Year); DOY <- c()

for(i in 1:length(Yrs)){
  tmpdates <- PhenoSite_0$Date[PhenoSite_0$Year==Yrs[i]]
  Jan1Jul <- as.Date(paste(Yrs[i],"-01-01", sep=""))
  ObsJulDayYr <- julian(as.Date(as.character(tmpdates),"%m/%d/%Y"), origin=Jan1Jul)
  DOY <- c(DOY,ObsJulDayYr)
}

PhenoSite_0 <- cbind(PhenoSite_0, DOY)
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

#Fits model per species using snow disappearance as prediction - across trails / plots / years
#haven't tried this yet - needs trouble shooting
curvefit_snowmod <- function (param){ #curve fitting function for mle
  intpeak  <- param[1]
  slopepeak <- param[2]
  rangep <- param[3]
  maxp   <- param[4]
  peakp <- intpeak + slopepeak*SDD
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


#################################################
######### Fit plot specific model  ##############
#################################################

#Define species of interest (for reference, RL and GB focal species listed)
species <- c("ANOC","ARLA","ASLE","CAMI","CAPA","ERGR","ERMO","ERPE","LUAR","PEBR","VASI"); nspp <- length(species)
#REFLECTION LAKES SPECIES ARE <- c("ANOC","CAPA","ERMO","ERPE","LIGR","LUAR","MIAL","PEBR","POBI","VASI")
#GLACIER BASIN SPECIES ARE <- c("ANAR","ARLA","ASLE","CAMI", "ERGR","LUAR","MEPA","PEBR","POBI","VASI")
#Exclude non focal species of interest
specieskeep <- PhenoSite$Species %in% species
PhenoSite_Focal <- PhenoSite[specieskeep,]

#how many years?
years <- unique(PhenoSite_Focal$Year)

#where to save output
parameters <- c() #save parameters and SDD from plot specific fits

##Now nested for loops to fit curves for all year, species, plots
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
      tmp_pars <- c(years[i], plots[j], SDDplt, speciesinplot[k], model1$par[1:3])
      parameters <- rbind(parameters, tmp_pars)
    }
  }
}
  
#turn parameters into a data frame
dimnames(parameters) <- list(c(), c("year","plot","SDD","species","peak","range","max"))
parameters <- data.frame(parameters)

#change storage type to numeric - all except plot (since a few plots have a, b)
parameters$year <- as.numeric(parameters$year)
parameters$species <- as.factor(parameters$species)
parameters$SDD <- as.numeric(parameters$SDD)
parameters$peak <- as.numeric(parameters$peak)
parameters$range <- as.numeric(parameters$range)
parameters$max <- as.numeric(parameters$max)

##############
####Various plots - by year, year / plot, by year / species###

#define plotting colors by species
plotcol <- c("pink","orangered","yellow","purple","lightblue","grey","magenta","yellowgreen",
             "navyblue","azure4","yellow4","yellowgreen","orchid","turquoise","salmon")

##FIRST PLOT - all curves, all years; per trail

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
    text(c(120,150,180,210,240,270),-0.1, labels=c("April","May","Jun","Jul","Aug","Sep"))
    text(110,0.9,labels=years[i])
  
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

##NEXT PLOT - plot curves separate by year (1 plot per year)
#plotting parameters
years <- unique(parameters$year)

#one plot per year
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
    text(c(150,180,210,240,270),-0.075, labels=c("May","June","July","August","September"))
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


##Third plot - curves per species & year (one plot per year)
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
    text(c(150,180,210,240,270),-0.075, labels=c("May","Jun","Jul","Aug","Sep"))
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

##Plot SDD vs. optim, range, per species
for(i in 1:length(species)){
  parspecies <- parameters[parameters$species==species[i],]
  X11(width=7.5, height=4.5)
  par(mfrow=c(1,2),omi=c(0,0,0,0), mai=c(0.35,0.335,0.3,0.0),tck=-0.02, mgp=c(1.25,0.25,0),xpd=TRUE)
  
  plot(parspecies$SDD,parspecies$peak, xlab="SDD",ylab="peak", pch=21, bg=plotcol[i])
  title(paste(species[i],"- SDD vs. optim"))
  
  plot(parspecies$SDD,parspecies$range, xlab="SDD",ylab="peak", pch=21, bg=plotcol[i])
  title(paste(species[i],"- SDD vs. range"))
  
}       
        

###Snippet of code that could be useful when fitting model per trail / year
#Pull out data for one of two transects - uncomment statement
PhenoSite_Trail <- PhenoSite[PhenoSite$Transect.x=="Reflection Lakes",]
#PhenoSite_Trail <- PhenoSite[PhenoSite$Transect.x=="Glacier Basin",]

#Define species of interest (for reference, RL and GB focal species listed)
species <- c("ANOC","CAPA","ERMO","ERPE","LUAR","PEBR","VASI"); nspp <- length(species)
#REFLECTION LAKES SPECIES ARE <- c("ANOC","CAPA","ERMO","ERPE","LIGR","LUAR","MIAL","PEBR","POBI","VASI")
#GLACIER BASIN SPECIES ARE <- c("ANAR","ARLA","ASLE","CAMI", "ERGR","LUAR","MEPA","PEBR","POBI","VASI")
#Exclude non focal species of interest
specieskeep <- PhenoSite_Trail$Species %in% species
PhenoSite_Trail <- PhenoSite_Trail[specieskeep,]

