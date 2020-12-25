########################################
## Analysis for MeadoWatch Data
## Goals: To determine the relationship between snow disappearance data and phenology
##        for the 16 focal MeadoWatch species (Reflection Lakes and Glacier Basin Transect)
## This script: 
##     A. Fits curve to each year / plot / species / phenophase combo
##     B. Creates plots of data + curves for all years / species / phenophase (to check)
##     C. Creates figures of general interest (for orientation, etc)
## Last worked on by Janneke Hille Ris Lambers: April 21, 2019
#######################################


###############################################################
#########Read in Data, assemble important explanatory variables
###############################################################

# Read in data: Phenodat is phenology data; stationdat is lat / long and includes snow disappearance info
PhenoDat <- read.csv("MW_PhenoDat_2013_2018.csv", header=TRUE) #phenology data
StationDat <- read.csv("MW_SiteDat_2013_2018.csv", header=TRUE) #information about station
MergePhenoStation <- merge(StationDat,PhenoDat, by=c("Year","Site_Code"))

#Tidy Data into columns needed; order by year
PhenoSite_0 <- MergePhenoStation[,c(1:2, 5:6, 8, 11, 14:15, 17:19, 21, 23, 25, 27)]
PhenoSite_0 <- PhenoSite_0[order(PhenoSite_0[,1]),]

#Calculate Julian Days of observations; DSS=days since snow; add to PhenoSite
Yrs <- unique(PhenoSite_0[,1]); DOY <- c()

for(i in 1:length(Yrs)){
  tmpdates <- PhenoSite_0[PhenoSite_0[,1]==Yrs[i],6]
  Jan1Jul <- as.Date(paste(Yrs[i],"-01-01", sep=""))
  ObsJulDayYr <- julian(as.Date(as.character(tmpdates),"%m/%d/%Y"), origin=Jan1Jul)
  DOY <- c(DOY,ObsJulDayYr)
}

DSS <- DOY - PhenoSite_0$SDD
PhenoSite_0 <- cbind(PhenoSite_0, DOY, DSS)
PhenoSite <- PhenoSite_0[,c(1:6,16:17,7:15)] #reorganize data

############################################
#### Extract information on # volunteers####
############################################

PhenoSite_vol <- PhenoSite[PhenoSite$QA.QC==0,]
yrs <- unique(PhenoSite$Year)
sites <- unique(substr(PhenoSite$Site_Code,1,2))
totalhikes <- 0; totalobs <-0

for(i in 1:length(yrs)){
  PhenoSite_volyr <- PhenoSite_vol[PhenoSite_vol$Year==yrs[i],]
  for(j in 1:2){
    yearsite <- PhenoSite_volyr[substr(PhenoSite_volyr$Site_Code,1,2)==sites[j],]
    if(dim(yearsite)[1]==0){next}
      Obs <- paste(yearsite$Date,yearsite$Observer)
      nhikes <- length(unique(Obs))
      nobs <- dim(yearsite)[1]
      print(paste(yrs[i],sites[j]))
      print("number hikes"); print(nhikes)
      print("observations"); print(nobs)
      totalhikes <- totalhikes + nhikes
      totalobs <- totalobs + nobs
  }
}

###################################################################
#### Determine species in each plot (for data sheets, pamphlets####
###################################################################

yr <- 2018 #relevant year
Ph_RL_yr <- PhenoSite[substr(PhenoSite$Site_Code,1,2)=="RL" & PhenoSite$Year==yr,]
Ph_RL_yr[is.na(Ph_RL_yr)==TRUE] <- 0
Ph_RL_yr <- droplevels(Ph_RL_yr)

Ph_GB_yr <- PhenoSite[substr(PhenoSite$Site_Code,1,2)=="GB" & PhenoSite$Year==yr,]
Ph_GB_yr[is.na(Ph_GB_yr)==TRUE] <- 0
Ph_GB_yr <- droplevels(Ph_GB_yr)


#Determine for RL - what has been observed and how many
Obs_RL <- tapply(Ph_RL_yr$Year, list(Ph_RL_yr$Site_Code, Ph_RL_yr$Species),length)
bud_RL <- tapply(Ph_RL_yr$Bud, list(Ph_RL_yr$Site_Code, Ph_RL_yr$Species),sum)
flwr_RL <- tapply(Ph_RL_yr$Flower, list(Ph_RL_yr$Site_Code, Ph_RL_yr$Species),sum)
fruit_RL <- tapply(Ph_RL_yr$Fruit, list(Ph_RL_yr$Site_Code, Ph_RL_yr$Species),sum)
seed_RL <- tapply(Ph_RL_yr$Disperse, list(Ph_RL_yr$Site_Code, Ph_RL_yr$Species),sum)
bud_RL[is.na(bud_RL)==TRUE] <- 0; fruit_RL[is.na(fruit_RL)==TRUE] <- 0
flwr_RL[is.na(flwr_RL)==TRUE] <- 0; seed_RL[is.na(seed_RL)==TRUE] <- 0

all_RL <- bud_RL + flwr_RL + fruit_RL + seed_RL
prop_RL <- 0.25* all_RL / Obs_RL

#Determine for GB - what has been observed and how many
Obs_GB <- tapply(Ph_GB_yr$Year, list(Ph_GB_yr$Site_Code, Ph_GB_yr$Species),length)
bud_GB <- tapply(Ph_GB_yr$Bud, list(Ph_GB_yr$Site_Code, Ph_GB_yr$Species),sum)
flwr_GB <- tapply(Ph_GB_yr$Flower, list(Ph_GB_yr$Site_Code, Ph_GB_yr$Species),sum)
fruit_GB <- tapply(Ph_GB_yr$Fruit, list(Ph_GB_yr$Site_Code, Ph_GB_yr$Species),sum)
seed_GB <- tapply(Ph_GB_yr$Disperse, list(Ph_GB_yr$Site_Code, Ph_GB_yr$Species),sum)
bud_GB[is.na(bud_GB)==TRUE] <- 0; fruit_GB[is.na(fruit_GB)==TRUE] <- 0
flwr_GB[is.na(flwr_GB)==TRUE] <- 0; seed_GB[is.na(seed_GB)==TRUE] <- 0

all_GB <- bud_GB + flwr_GB + fruit_GB + seed_GB
prop_GB <- 0.25* all_GB / Obs_GB

###############################################################
#########Define Maximum Likelihood Models (to fit to data)
###############################################################

#Null model - just says that probability of flowering does not vary with time
nullfit <- function (param){ 
	meanp  <- param[1]
	pred   <- rep(meanp, times=length(phenophase))
	llik     <- dbinom(phenophase,1,pred, log=TRUE)
	return(-sum(llik))
}

#Alternative model - just says that probability of flowering varies with days (could be DSS, or day of year)
curvefit <- function (param){ #curve fitting function for mle
	 peakp  <- param[1]
	 rangep <- param[2]
	 maxp   <- param[3]
	 pred   <- maxp*(1/(rangep*sqrt(2*pi)))*exp(-0.5*((days-peakp)/rangep)^2)
	 llik     <- dbinom(phenophase,1,pred, log=TRUE)
	 return(-sum(llik))
}

#Predicts probability of observing phenological phase, given parameters (used for plotting)
predphen <- function (xx, param){ 
  days    <- xx
	peakp  <- param[1]
	rangep <- param[2]
	maxp   <- param[3]
	pred   <- maxp*(1/(rangep*sqrt(2*pi)))*exp(-0.5*((days-peakp)/rangep)^2)
	return(pred)
}


#######################################################
#########Fit models to data, save output ##############
#######################################################

#First define parameters
species <- unique(PhenoSite[,11]) #first identify unique species
species <- species[order(species)]
years <- unique(PhenoSite[,1])
sites <- unique(PhenoSite[,2])

#colors for plotting
plotcol <- c("azure3","yellow4","orchid1","orange","magenta","yellow","yellowgreen","pink",
             "gray","azure4","plum","cadetblue1","gold","darkblue","maroon","deeppink","lightblue")

#where to save output
output1 <- c() #timing of phenophase; parameter for model fit
output2 <-c() #model fit

##set plotting, DSS vs. DOY for time; phenophase
plotall <- 0 #set for no plotting; 1 if you want plots of raw data and fits (in phenophase directories, by year) to be saved
plotsnow <- 0 # 1 for days since snow, 0 for doy
phenocats <- c(13,14,15,16) # 13=buds; 14=flwrs; 15=fruit; 16=seeds
names(phenocats) <- c("bud","flwr","fruit","seed")

##Fit phenological curves to each species / site / year combo; with a for loop
for(i in 1:length(species)){ #loop for each species
  #extract data for that species
  speciesdat<-PhenoSite[PhenoSite[,11]==species[i],]
  
  #loop for each year
  for(j in 1:length(years)){ # loop for each year
    speciesyeardat <- speciesdat[speciesdat[,1]==years[j],]
    mxplots <- length(unique(speciesyeardat$Site_Code))
    rwcl <- ceiling(sqrt(mxplots)) 

    #break to next if no data for that species / year
    if(dim(speciesyeardat)[1]==0){next}
    
    for(k in 1:length(phenocats)){ # loop for each site
      #one plot for each site / year; if plotall is 1
      if(plotall==1){
        localdir <- getwd()
        if(dir.exists(paste(localdir,"/",names(phenocats)[k],sep=""))==FALSE)
           {dir.create(paste(localdir,"/",names(phenocats)[k],sep=""))}
        
      jpeg(file=paste(names(phenocats)[k],"/",species[i],"-",years[j],".jpg",sep=""), width=8,height=8,units="in", res=300)
      par(mfrow=c(rwcl,rwcl), omi=c(0,0,0,0), mai=c(0.25,0.25,0.25,0.25), mgp=c(0.75,0.1,0), tck=-0.01)}
      
      for(m in 1:length(sites)){ #loop for each phenological stage
        #extract data from species dat
        speciesyearsitedat <- speciesyeardat[speciesyeardat[,2]==sites[m],]
        if(dim(speciesyearsitedat)[1]==0){next}
      
        #define parameters for model fitting: days / time
        if(plotsnow==1){days <- speciesyearsitedat[,8]} #explanatory variable: DOY = column 7, DSS = column 8
        if(plotsnow==0){days <- speciesyearsitedat[,7]} #explanatory variable: DOY = column 7, DSS = column 8
      
        #define phenophase
        phenophase <- speciesyearsitedat[,phenocats[k]] #yes / no for phenophase    
      
        #remove days when no observations were made; NA in phenophase
	      days <- days[is.na(phenophase)==FALSE]
	      phenophase <- phenophase[is.na(phenophase)==FALSE]
	    
	      #go to next if species was not observed more than 2 times during the summer
	      if(length(phenophase[phenophase[]==1])<3){next}
	    
	      #add three weeks of zeroes (1 week, 2 weeks, and 3 weeks prior to SDD)
	      mxdays <- max(days)
	      SDDsiteyear <- speciesyearsitedat[1,5]
	      if(plotsnow==1){
	        days <- c(-21, -14, -7,days,mxdays+7, mxdays+14, mxdays+21)
	      }
	      if(plotsnow==0){
	        days <- c(SDDsiteyear-21, SDDsiteyear-14, SDDsiteyear-7,days,mxdays+7, mxdays+14, mxdays+21)
	      }
	      phenophase <- c(0,0,0,phenophase,0,0,0)

        #now fit model, save data
        model0 <- optimize(nullfit, c(0.000001,0.999999)) #fit null model
	
      	#now fit curve
      	param <- c(mean(days), sd(days), 0.5) # initial parameters: model fits pretty fussy about this
      	model1 <- optim(param, curvefit, control = list(maxit = 20000))
     	  if(model1$convergence==1){print(paste(species[i],years[j],sites[m], sep="-"))}
    	
    	  #calculate AIC, p value
    	  AICnull <- round(2*(model0$objective+1),1)
    	  AICalt <- round(2*(model1$value + 3),1)
    	  pcurve <- signif(pchisq(model1$value-model0$objective,2),3)
	
    	  #add to output tables: parameters
    	  tmpoutput1 <- c(as.character(species[i]), years[j], as.character(sites[m]), 
    	                  as.character(names(phenocats)[k]), SDDsiteyear, round(model1$par, 4))
    	  output1 <- rbind(output1, tmpoutput1)
    	
    	  #add output to output tables: model fitting
      	tmpoutput2 <- c(as.character(species[i]), years[j], as.character(sites[m]), 
    	                  as.character(names(phenocats)[k]), AICnull, AICalt, pcurve)
      	output2 <- rbind(output2, tmpoutput2)
  
    	  #calculate predicted based on model3
       	xx <- seq(min(days), max(days))
       	yy <- predphen(xx, model1$par)

     	  #now plot; if plotall is 1: need to rethink plotting
     	  if(plotall==1){
         	ifelse(plotsnow==1,xlbl<-"Days since Snow", xlbl<-"Julian Days")
        	plot(xx, yy, xlab=xlbl, ylab="p(pheno)",xlim=c(min(days), max(days)), ylim=c(0,1),
      	       type="l", col=plotcol[i], lwd=2)
         	points(days, phenophase, pch=21, bg=plotcol[i])
     	    title(paste(species[i],"-",years[j],names(phenocats)[k],"-",sites[m], sep=""), cex=0.5)
     	  }
      }
      if(plotall==1){dev.off()} #turns the device 'off, necesssary for plotting
    }
  }
}

#Output tables - add dimension names, make appropriate columns numeric
dimnames(output1) <- list(c(), c("Species","Year","Site","Phenophase","SDD", "Opt","Duration","Height"))
outputpars <- data.frame(output1)
for(i in 5:8){outputpars[,i] <- as.numeric(as.character(outputpars[,i]))}

#Output tables - add dimension names, make appropriate columns numeric
dimnames(output2) <- list(c(), c("Species","Year","Site","Phenophase","AICnull", "AICalt","pcurve"))
outputfit <- data.frame(output2)
for(i in 5:7){outputfit[,i] <- as.numeric(as.character(outputfit[,i]))}


#########################################################
############# Plot SDD vs. opt, all phenophases ########
#########################################################

## To do: automate adding of color by year

#Analyze and plot DSS vs. optimum: all 4 phenophases (on one graph; one graph per species)
for(i in 1:length(species)){
  sppdat <- outputpars[outputpars[,1]==species[i],]

  #One plot per species - separate for peak and duration
  print(species[i])
  X11(width=6,height=6)
  par(mfrow=c(2,2), omi=c(0,0,0,0), mai=c(0.4,0.4,0.4,0.4),tck=-0.02, mgp=c(1.25,0.5,0))
  
  for(j in 1:length(phenocats)){ #now by phenophase
    sppphenodat <- sppdat[sppdat[,4]==names(phenocats)[j],]
    
    #go to next if too few observations; less than 2
    if(dim(sppphenodat)[1]<2){next}
    
    #define parameters
    DSS <- sppphenodat[,5]; peakp <- sppphenodat[,6]
    
    #define color, shape for plotting symbols (shape = transect, color = year)
    plotshp <- rep(21, length=dim(sppphenodat)[1])
    siteloc <- substring(sppphenodat[,3],1,2)
    plotshp[siteloc[]=="GB"] <- 24
    plotcol <- rep("yellowgreen", length=dim(sppphenodat)[1])
    plotcol[sppphenodat[,2]==2014] <- "magenta"
    plotcol[sppphenodat[,2]==2015] <- "orange"
    plotcol[sppphenodat[,2]==2016] <- "navyblue"
    plotcol[sppphenodat[,2]==2017] <- "yellow"
    plotcol[sppphenodat[,2]==2018] <- "purple"
    
    #Analysis: is peak flowering a function of DSS?
    DSStest <- lm(peakp~DSS)

    #Make a plot
    plot(DSS, peakp, xlab="Snowmelt (Julian Days)", ylab=paste("peak(",names(phenocats)[j],")",sep=""),
         pch=plotshp, bg=plotcol, cex=1.5)
    abline(0,1)
    abline(coef(DSStest)[1], coef(DSStest)[2], lty=2)
    title(paste(species[i], names(phenocats)[j], sep="-"))
    if(j==1){legend(x="topleft", legend=c("GB","RL","2013","2014","2015","2016","2017","2018"), pch=c(24,21,21,21,21,21,21,21), 
                    pt.bg=c("grey","grey","yellowgreen","magenta","orange","navyblue","yellow","purple"), cex=0.8)}
  }
}


#####################################################################
############# Plot DSS vs opt, DSS vs. duration: flowering ########
##################################################################

DSSoutput <- array(data=NA, dim=c(length(species), 4, 2), 
                   dimnames=list(species,  c("mean","intercept","slope","p-value"), c("optim","duration")))

#Analyze and plot DSS vs. optimum: all 4 phenophases (on one graph; one graph per species)
for(i in 1:length(species)){
  sppdat <- outputpars[outputpars[,1]==species[i],]
  
  #One plot per species - separate for peak and duration
  print(species[i])
  X11(width=9,height=5)
  par(mfrow=c(1,2), omi=c(0,0,0,0), mai=c(0.4,0.4,0.4,0.4),tck=-0.02, mgp=c(1.25,0.5,0))
  sppphenodat <- sppdat[sppdat[,4]=="flwr",]
    
  #go to next if too few observations; less than 2
  if(dim(sppphenodat)[1]<2){next}
    
  #define parameters
  DSS <- sppphenodat[,5]; peakp <- sppphenodat[,6]; durp <- sppphenodat[,7]
    
  #define color, shape for plotting symbols (shape = transect, color = year)
  plotshp <- rep(21, length=dim(sppphenodat)[1])
  siteloc <- substring(sppphenodat[,3],1,2)
  plotshp[siteloc[]=="GB"] <- 24
  plotcol <- rep("yellowgreen", length=dim(sppphenodat)[1])
  plotcol[sppphenodat[,2]==2014] <- "magenta"
  plotcol[sppphenodat[,2]==2015] <- "orange"
  plotcol[sppphenodat[,2]==2016] <- "navyblue"
  plotcol[sppphenodat[,2]==2017] <- "yellow"
  plotcol[sppphenodat[,2]==2018] <- "purple"
    
  #Analysis: is peak flowering a function of DSS?
  nulltestpeak <- lm(peakp~1)
  DSStestpeak <- lm(peakp~DSS)
  DSSoutput[i,,1]<-c(nulltestpeak$coefficients, DSStestpeak$coefficients, summary(DSStestpeak)[[4]][2,4])
    
  #Analysis: is flowering duration a function of DSS?
  nulltestdur <- lm(durp~1)
  DSStestdur <- lm(durp~DSS)
  DSSoutput[i,,2]<-c(nulltestdur$coefficients, DSStestdur$coefficients, summary(DSStestdur)[[4]][2,4])
    
  #Make a plot for peak flowering as a function of DSS
  plot(DSS, peakp, xlab="Snowmelt (Julian Days)", ylab="peak flowering",
       pch=plotshp, bg=plotcol, cex=1.5)
  abline(0,1)
  abline(coef(DSStestpeak)[1], coef(DSStestpeak)[2], lty=2)
  title(paste(species[i], "peak flowering", sep="-"))
  legend(x="topleft", legend=c("GB","RL","2013","2014","2015","2016","2017","2018"), pch=c(24,21,21,21,21,21,21,21), 
                    pt.bg=c("grey","grey","yellowgreen","magenta","orange","navyblue","yellow","purple"), cex=0.8)
  mtext(paste("p(snow)=",round(summary(DSStestpeak)[[4]][2,4],3)), side=3, line=-1, adj=0.5)

  
  #Make a plot for duration flowering as a function of DSS
  plot(DSS, durp, xlab="Snowmelt (Julian Days)", ylab="duration flowering",
       pch=plotshp, bg=plotcol, cex=1.5)
  abline(coef(DSStestdur)[1], coef(DSStestdur)[2], lty=2)
  title(paste(species[i], "duration flowering", sep="-"))
  mtext(paste("p(snow)=",round(summary(DSStestdur)[[4]][2,4],3)), side=3, line=-1, adj=0.5)
  
}

###Now create a plot that would help identify richness, duration of different years
DSSobs <- StationDat$SDD[StationDat$Transect=="Reflection Lakes"]
minDSS <- min(DSSobs); maxDSS <- max(DSSobs)+5 #x = 120-210; y=150-250
X11(width=9,height=9)
#flwrcol <- c("grey","yellow","papayawhip","pink","orange","lightblue","goldenrod","magenta","brown",
#             "purple","cornsilk2","yellowgreen","salmon","sienna","rosybrown", "turquoise","violetred")
species_RL <- c("ANOC","CAPA","ERMO","ERPE","LIGR","LUAR","PEBR","POBI","VASI")
flwrcolRL <- c("yellowgreen","magenta","yellow","pink","grey","purple","lightblue","maroon","seagreen")

for(i in 1:length(species_RL)){
  if(i==1){
    plot(c(150, 230), xlim=c(minDSS, maxDSS), ylim=c(150,260), type="n", xlab="Days since snow", ylab="Peak Flowering") 
  
  #also plot lines indicating range of SDD for years
  yrs <- seq(2013,2018)
  for(j in 1:length(yrs)){
    yrDSS <- StationDat[StationDat$Year==yrs[j]&StationDat$Transect=="Reflection Lakes",]
    minDSSyr <- min(yrDSS$SDD); maxDSSyr <- max(yrDSS$SDD)
    linepos <- 260 + (2013 - yrs[j])*3
    lines(c(minDSSyr,maxDSSyr), c(linepos,linepos), lwd=2, col="lightblue")
    text((maxDSSyr+5),linepos, labels=yrs[j])
    #polygon(c(minDSSyr,minDSSyr,maxDSSyr,maxDSSyr), c(150,250,250,150),density=10,col="lightblue")
  }}
  
  
  DSSoutput_spp <- DSSoutput[which(dimnames(DSSoutput)[[1]]==species_RL[i]),,]
  int <- DSSoutput_spp[2,1]; slope <- DSSoutput_spp[3,1]
  abline(int,slope, col=flwrcolRL[i], lwd=2)

}

###Now create a plot that would help identify duration per plot
DSSobs <- StationDat$SDD[StationDat$Transect=="Reflection Lakes"]
minDSS <- min(DSSobs); maxDSS <- max(DSSobs)+5 #x = 120-210; y=150-250
X11(width=9,height=9)
#flwrcol <- c("grey","yellow","papayawhip","pink","orange","lightblue","goldenrod","magenta","brown",
#             "purple","cornsilk2","yellowgreen","salmon","sienna","rosybrown", "turquoise","violetred")
species_RL <- c("ANOC","CAPA","ERMO","ERPE","LIGR","LUAR","PEBR","POBI","VASI")
flwrcolRL <- c("yellowgreen","magenta","yellow","pink","grey","purple","lightblue","maroon","seagreen")

for(i in 1:length(species_RL)){
  if(i==1){
    plot(c(150, 230), xlim=c(minDSS, maxDSS), ylim=c(150,260), type="n", xlab="Days since snow", ylab="Peak Flowering") 
    
    #also plot points indicating snow duration in 2016
    yrDSS <- StationDat[StationDat$Year==2016&StationDat$Transect=="Reflection Lakes",]
    points(yrDSS$SDD,rep(250,length=dim(yrDSS)[1]), pch=21, cex=2,bg="lightblue")
    text(yrDSS$SDD,rep(255,length=dim(yrDSS)[1]), labels=yrDSS$Site_Code)
    }
  
  
  DSSoutput_spp <- DSSoutput[which(dimnames(DSSoutput)[[1]]==species_RL[i]),,]
  int <- DSSoutput_spp[2,1]; slope <- DSSoutput_spp[3,1]
  abline(int,slope, col=flwrcolRL[i], lwd=2)
  
}





###############################################################
############## Figures for various publications ###############
##############################################################


####################################################
##Plot date vs flowering and seeding - LUAR (VASI for Mountain views)
sppdat <- outputpars[outputpars$Species=="LUAR",]
sppdatflwr <- sppdat[sppdat$Phenophase=="flwr",]
sppdatsd <- sppdat[sppdat$Phenophase=="seed",]
DSSflwr <- sppdatflwr[,5]; peakpflwr <- sppdatflwr[,6]
DSSsd <- sppdatsd[,5]; peakpsd <- sppdatsd[,6]

#define color, shape for plotting symbols (shape = transect, color = year)
plotshpflwr <- rep(21, length=DSSflwr)
plotshpsd <- rep(21, length=DSSsd)
plotcolflwr <- rep("blueviolet", length=DSSflwr)
#plotcolflwr[sppdatflwr[,2]==2015] <- "orange"
plotcolsd <- rep("darkgoldenrod4", length=DSSsd)
#plotcolsd[sppdatsd[,2]==2015] <- "orange"

X11(width=5,height=5)
par(mfrow=c(1,1), omi=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.5),tck=-0.02, mgp=c(1.15,0.35,0),xpd=TRUE)

#Make a plot
plot(DSSflwr, peakpflwr, xlab="Snowmelt (Julian Days)", ylab="Peak Flowering or Seeding",
     pch=plotshpflwr, bg=plotcolflwr, cex=1.5, xaxt="n", yaxt="n",
     xlim=c(100,200), ylim=c(165,270))
points(DSSsd,peakpsd, pch=plotshpsd, bg=plotcolsd, cex=1.5)

#add bestfit lines
#lmflwr <- lm(peakpflwr~DSSflwr); abline(coef(lmflwr), col="blueviolet")
#lmsd <- lm(peakpsd~DSSsd); abline(coef(lmsd), col="darkgoldenrod4")

lines(c(162,203),c(162,203), col="lightblue", lwd=2, lty=1)
text(c(105,135,165,195),159,labels=c("April","May","June","July"))
text(94, c(165,195,225,255),labels=c("June","July","Aug","Sept"), srt=90)

title("Subalpine lupine: 2013-2018")
legend(x="topleft", legend=c("flower","seed"), pch=c(21,21), pt.cex=1.5,
                pt.bg=c("blueviolet", "darkgoldenrod4"))


#################################################################
##### Plot Flowering Season, average per year 
yrs <- unique(outputpars$Year); yrs <- yrs[order(yrs)]
sites <- c("GB","RL")
phen <- "flwr" #which phenophase
pltcls <- c("limegreen","yellowgreen","yellow","lightpink","orange","magenta","goldenrod","yellow2",
            "pink","grey","grey","purple","royalblue","yellow3","springgreen","maroon","lightblue")

for(i in 1:length(yrs)){
  outputyr <- outputpars[outputpars$Year==yrs[i],]
  outputyr <- outputyr[outputyr$Phenophase==phen,]
  
  #set plotting functions
  X11(width=9,height=4)
  par(mfrow=c(1,2), mgp=c(1.2,0.5,0), tck=-0.02)
  
  for(j in 1:2){
    outputsite <- outputyr[substr(outputyr$Site,1,2)==sites[j],]
    meanopt <- tapply(outputsite$Opt, outputsite$Species, mean)
    meandur <- tapply(outputsite$Duration, outputsite$Species, mean)
    meanhgt <- tapply(outputsite$Height, outputsite$Species, mean)
    
    #Now set up plot
    time <- seq(150,270)
    plot(time,rep(1,times=length(time)),type="n",ylim=c(0,2),xlab="Time",ylab="P(flower)",xaxt="n")
    
    for(k in 1:length(meanopt)){
      if(names(meanopt)=="LICA"){next}
      if(is.na(meanopt[k])==TRUE){next}
      prs <- c(meanopt[k],meandur[k],meanhgt[k])
      flr <- predphen(time,prs)
      time2 <- time[flr[]>0.01]
      flr2 <- flr[flr[]>0.01]
      lines(time2,flr2,col=pltcls[k])
    }
    title(paste(yrs[i],sites[j]))
  }
}


#################################################################
#####Flowering Richness  vs time (empirical)
ObsDat <- paste(PhenoSite$Year, substr(PhenoSite$Site_Code,1,2), PhenoSite$DOY, PhenoSite$Observer) 
ObsHikes <- unique(ObsDat) #unique observers by year and date
HikingExperience <- c()

for(i in 1:length(ObsHikes)){
  HikeDat <- PhenoSite[ObsDat[]==ObsHikes[i],]
  hikeinfo <- HikeDat[1,1:9]
  HikeDat <- HikeDat[is.na(HikeDat$Flower)==FALSE,]
  NSppFlwr <- length(unique(HikeDat$Species[HikeDat$Flower==1]))
  NPltsFlwr <- length(unique(HikeDat$Site_Code[HikeDat$Flower==1]))
  plotbyf <- tapply(HikeDat$Flower, list(droplevels(HikeDat$Site_Code),droplevels(HikeDat$Species)), sum)
  totplotf <- rowSums(plotbyf, na.rm=TRUE)
  AvgFlwrs <- mean(totplotf); if(length(totplotf)==0){Avgflwrs<-0}
  MaxFlwrs <- max(totplotf); if(length(totplotf)==0){Maxflwrs<-0}
  #Now assemble data frame

  hikeinfo$NSppFlwr <- NSppFlwr
  hikeinfo$NPltsFlwr <- NPltsFlwr
  hikeinfo$AvgFlwrs <- AvgFlwrs
  hikeinfo$MaxFlwrs <- MaxFlwrs
  HikingExperience <- rbind(HikingExperience, hikeinfo)
}

#Change site_code to site
HikingExperience$Site_Code <- substr(HikingExperience$Site_Code,1,2)
dimnames(HikingExperience)[[2]][2] <- "Site"

#now plot, by year and trail
yrs <- unique(HikingExperience$Year); yrs <- yrs[order(yrs)]
sites <- c("GB","RL")
season <- c(min(HikingExperience$DOY),max(HikingExperience$DOY))

for(i in 1:length(yrs)){
  outputyr <- HikingExperience[HikingExperience$Year==yrs[i],]

  #set plotting functions
  X11(width=8,height=9)
  par(mfcol=c(4,2), mgp=c(1.2,0.5,0), tck=-0.02, omi=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.3))

  #Now add lines
    
  for(j in 1:2){
    outputsite <- outputyr[outputyr$Site==sites[j],]
    if(dim(outputsite)[1]==0){next}
   
    #plot total flowers, plots with flowers, average flowers per plot, max flowers per plot
    plot(outputsite$DOY,outputsite$NSppFlwr, pch=21, bg="pink",
         xlim=season, ylim=c(0,11), xlab="Time",ylab="Species Flowering")
    title(paste(yrs[i], sites[j],"Total Species Flowering"))
    
    plot(outputsite$DOY,outputsite$NPltsFlwr,pch=21, bg="pink",
         xlim=season, ylim=c(0,18),xlab="Time",ylab="Plots with Flowers")
    title(paste(yrs[i], sites[j],"Plots with Flowers"))

    plot(outputsite$DOY,outputsite$AvgFlwrs,pch=21, bg="pink",
         xlim=season, ylim=c(0,8), xlab="Time",ylab="Flowering Spp per Plot")
    title(paste(yrs[i], sites[j],"Species per Plot"))
    
    plot(outputsite$DOY,outputsite$MaxFlwrs,pch=21, bg="pink",
         xlim=season, ylim=c(0,10), xlab="Time",ylab="Max Flowers per Plot")
    title(paste(yrs[i], sites[j],"Max Species Flowering"))
    
  }
}

###################################
##Pretty plot, lines for 2018 newsletter
X11(width=5,height=6)
par(mfcol=c(2,1), mgp=c(1.0,0.25,0), tck=-0.02, omi=c(0,0,0,0), mai=c(0.5,0.5,0.4,0.3),xpd="TRUE")
yearcols <- c("forestgreen","royalblue","coral","pink","darkorchid1","cyan")
sitenames <- c("Glacier Basin", "Reflection Lakes")

for(i in 1:2){
  outputsite <- HikingExperience[HikingExperience$Site==sites[i],]
    plot(outputsite$DOY, outputsite$AvgFlwrs, type="n",
         xlim=c(160,260), ylim=c(0,4), xaxt="n", xlab="Date",ylab="Flowering Species per Plot")
    text(c(165,195,225,255),-0.5,labels=c("May","June","July","Aug"))
    title(sitenames[i])

    if(i==1){
      legend(x="topleft",legend=yrs[3:6],col=yearcols[3:6],cex=0.8, y.intersp=0.75, lwd=2)
    }
        
  for(j in 3:length(yrs)){
    outputyr <- outputsite[outputsite$Year==yrs[j],]
    if(dim(outputyr)[1]==0){next}

    #plot lines showing total species in flower
    #tmpxy <- cbind(outputyr$DOY,outputyr$AvgFlwrs)
    #tmpxy <- tmpxy[order(tmpxy[,1]),]
    #tmpxy2 <- cbind(tmpxy[,1],runmean(tmpxy[,2],20))
    #tmpxy2 <- tmpxy2[tmpxy2[,1]>160,]
    #tmpxy2 <- tmpxy2[tmpxy2[,1]<260,]
    #lines(tmpxy2[,1],tmpxy2[,2], col=yearcols[j], lwd=2)
    
    #plot lines showing total species in flower
    ltest <- loess(AvgFlwrs~DOY,outputyr, span=0.75)
    lpred <- predict(ltest, data.frame(DOY = seq(season[1],season[2],1)))
    lpred2 <- cbind(seq(season[1],season[2],1), lpred)
    lpred2 <- lpred2[is.na(lpred2[,2])==FALSE,]
    lpred2 <- lpred2[lpred2[,2]>0.01,]
    lpred2 <- lpred2[lpred2[,1]>160,]
    lpred2 <- lpred2[lpred2[,1]<260,]
    lines(lpred2, col=yearcols[j], lwd=3)
  }
}


#########################################################
############# Determine herbivory vs. flowering #
##First extract data from years where herbivory was sampled

PhenoSiteHerb <- PhenoSite[PhenoSite$Year>2016,] #2017 onwards
#PhenoSiteHerb <- PhenoSiteHerb[PhenoSiteHerb$Species=="ANOC"|PhenoSiteHerb$Species=="VASI"|PhenoSiteHerb$Species=="PEBR",]

#Determine unique year, site, plot, species combo
ObsPlot <- paste(PhenoSiteHerb$Year, PhenoSiteHerb$Site_Code,PhenoSiteHerb$Species)
TotObs <- unique(ObsPlot)

#specify output plot
herbdat <- c()

for(i in 1:length(TotObs)){
  outputplot <- PhenoSiteHerb[ObsPlot==TotObs[i],]
  if(dim(outputplot)[1]<5){next}
  pherb <- length(na.omit(outputplot$Herb))/length(outputplot$Flower)
  pflwr <- sum(na.omit(outputplot$Flower)) /length(outputplot$Flower) 
  tmp <- c(as.character(outputplot$Species[1]),substr(outputplot$Site_Code[1],1,2),pherb,pflwr)
  herbdat <- rbind(herbdat,tmp)
}

herbdat2 <- data.frame(herbdat)
dimnames(herbdat2)[[2]] <- c("spp","transect","herb","flwrobs")
herbdat2$flwrobs <- as.numeric(as.character(herbdat2$flwrobs))
herbdat2$herb <- as.numeric(as.character(herbdat2$herb))
summaryherb <- tapply(herbdat2$herb, list(herbdat2$spp,herbdat2$transect),mean)
summaryherb2 <- summaryherb[is.na(summaryherb[,1])==FALSE,]
summaryherb2 <- summaryherb2[is.na(summaryherb2[,2])==FALSE,]
summaryherb2 <- summaryherb2[-1,]


x11(width=5,height=3)
par(mgp=c(1.25,0.4,0), mai=c(0.4,0.5,0.4,0.3),omi=c(0,0,0,0), tck=-0.02)
bcol <- c("orchid3","yellow","maroon","yellowgreen")
barplot(summaryherb2,col=bcol,beside=TRUE, ylim=c(0,0.15), ylab="P(Floral Herbivory)", names.arg=c("Glacier Basin","Reflection Lakes") )
abline(h=0)
legend(x="topright", legend=c("Subalpine lupine","Bracted lousewort","American bistort","Sitka valerian"), cex=0.75, fill=bcol)


##############################
####Plot snowmelt vs. peak flowering: different species
###RL Species: VASI, POBI, PEBR, LIGR, ERPE, ERMO (?), CAPA, LUAR
###GB Species: MEPA, ERGR, CAMI, ARLA, ASLE

#spp <- c("ERMO","CAPA","ERPE","LIGR","LUAR","POBI","PEBR","VASI")
#sppcol <- c("goldenrod","magenta","pink","lightblue","purple","maroon",
#            "yellowgreen", "grey")
spp <- c("ERGR","LUAR","CAMI","ARLA","ASLE","MEPA","POBI","PEBR")
sppname <- c("Glacier lily","Subalpine lupine","Scarlet paintbrush",
             "Broadlead arnica", "Cascade aster","Tall bluebell",
             "American bistort","Bracted lousewort")
sppcol <- c("gold","purple","orangered","orange","pink","skyblue",
            "maroon","yellowgreen")
rngssnow <- c(min(outputpars$SDD[substr(outputpars$Site,1,2)=="GB"]), 
              max(outputpars$SDD[substr(outputpars$Site,1,2)=="GB"]))
rngsflwr <- c(min(outputpars$Opt[substr(outputpars$Site,1,2)=="GB"]), 
              max(outputpars$Opt[substr(outputpars$Site,1,2)=="GB"]))

X11(width=6, height=5)
par(omi=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.5), mgp=c(1.25,0.5,0), tck=-0.02, xpd=NA)

#plot, but no lines
plot(165,195, xlab="Date of Snowmelt",ylab="Peak Flowering",
     xlim=c(100,200), ylim=c(135,270), type="n", xaxt="n", yaxt="n", xpd=NA)
text(c(105,135,165,195),125,labels=c("April","May","June","July"))
text(92, c(135,165,195,225,255),labels=c("May","June","July","Aug","Sept"), srt=90)


for(i in 1:length(spp)){
  sppdat <- outputpars[outputpars$Species==spp[i]&outputpars$Phenophase=="flwr",]
  spptest <- lm(sppdat$Opt~sppdat$SD)
  xx <- c(100,200)
  yy <- c(spptest$coef[1]+spptest$coef[2]*xx[1], spptest$coef[1]+spptest$coef[2]*xx[2])
  lines(xx,yy, col=sppcol[i], lwd=2)
}
legend(x="topleft", legend=sppname, col=sppcol, lwd=2, cex=0.65)

lines(c(135,200),c(135,200), lwd=2, col="navyblue")

######################
##Grab bag of plots from previous presentations: not in very good shape
###########

#1. Individual species: observations vs. curve
sppplot <- "CAPA" #pick species to plot; check 'species' for possibilities
phenoplot <- "flwr" #pick phenophase; check 'phenocats'
yearplot <- 2018
plotplot <- "RL5"

#extract raw data
SppDat <- PhenoSite[PhenoSite$Year==yearplot & PhenoSite$Species==sppplot & PhenoSite$Site_Code == plotplot, ]
#First fit relationship, for all phenophases
outputspp <- c()

for(i in 13:16){
  days <- SppDat$DOY
  phenophase <- SppDat[,i]
  days <- days[is.na(phenophase)==FALSE]; phenophase <- phenophase[is.na(phenophase)==FALSE]
  param <- c(mean(days),10,0.9)
  model1 <- optim(param, curvefit, control = list(maxit = 20000))
  outputspp <- rbind(outputspp,model1$par)
}
dimnames(outputspp) <- list(c("bud","flower","fruit","seed"),c("peak","range","max"))


#Flower only
#png(filename="PEBRpnts.png",width=6,height=5, units="in", res=1200)
X11(width=5,height=5)
par(mfrow=c(1,1), omi=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.5), mgp=c(1.1,0.3,0), tck=-0.02, xpd="TRUE")
plot(SppDat$DOY, SppDat$Flower, pch=21, bg="darkmagenta", cex=1.5, ylim=c(0,1), xlim=c(180,270), 
     xlab="", ylab="Flower Observations", xaxp=c(180,240,2), xaxt="n", yaxt="n")
text(c(195,225,255),-0.075,labels=c("June","July","August"),cex=1.5)
text(172.5,c(0.025,0.975),labels=c("No","Yes"), srt=90,cex=1.5)
#dev.off()

#Flower and line
time <- seq(180,270, by=0.25)
mxsppltyr <- outputspp[2,3]; rangesppltyr <- outputspp[2,2]; peaksppltyr <- outputspp[2,1]
probf <- mxsppltyr*(1/(rangesppltyr*sqrt(2*pi)))*exp(-0.5*((time-peaksppltyr)/rangesppltyr)^2)

#
lines(time,probf, col="darkmagenta", lwd=2)
lines(c(peaksppltyr-0.5*rangesppltyr,peaksppltyr+0.5*rangesppltyr),c(0.5,0.5))

#Now plot
X11(width=5,height=5)
#png(filename="ERGRpntslns.png",width=6,height=5, units="in", res=600)
par(mfrow=c(1,1), omi=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.5), mgp=c(1.1,0.3,0), tck=-0.02, xpd="NA")
plot(time,probf, type="l", col="darkmagenta", lwd=2, xlab="", xaxt="n", yaxt="n", 
     ylab="Flower Observations", ylim=c(0,1), xaxp=c(150,240,3))
points(SppDat$DOY[SppDat$DOY[]<240], SppDat$Flower[SppDat$DOY[]<240], pch=21, bg="gold1",cex=1.5)
text(c(195,225,255),-0.075,labels=c("July","August","September"),cex=1.5)
text(172.5,c(0.025,0.975),labels=c("No","Yes"), srt=90,cex=1.5)
#dev.off()

#Next plot all phenophases
time <- seq(180,270, by=0.25)
probb <- outputspp[1,3]*(1/(outputspp[1,2]*sqrt(2*pi)))*exp(-0.5*((time-outputspp[1,1])/outputspp[1,2])^2)
probf <- outputspp[2,3]*(1/(outputspp[2,2]*sqrt(2*pi)))*exp(-0.5*((time-outputspp[2,1])/outputspp[2,2])^2)
probr <- outputspp[3,3]*(1/(outputspp[3,2]*sqrt(2*pi)))*exp(-0.5*((time-outputspp[3,1])/outputspp[3,2])^2)
probs <- outputspp[4,3]*(1/(outputspp[4,2]*sqrt(2*pi)))*exp(-0.5*((time-outputspp[4,1])/outputspp[4,2])^2)

X11(width=6,height=6)
#png(filename="ERGRphenos.png",width=6,height=5, units="in", res=1200)
par(mfrow=c(1,1), omi=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.5), mgp=c(1.1,0.3,0), tck=-0.02, xpd="NA")
plot(time[probb[]>0.005],probb[probb[]>0.005], type="l", col="yellowgreen", lwd=2, xlab="", xaxt="n",  
     ylab="P(Phenophase)", ylim=c(0,1), xlim=c(180,270),xaxp=c(150,270,4))
lines(time[probf[]>0.005],probf[probf[]>0.005], col="gold1", lwd=2)
lines(time[probr[]>0.005],probr[probr[]>0.005], col="green4", lwd=2)
lines(time[probs[]>0.005],probs[probs[]>0.005], col="darkgoldenrod", lwd=2)
text(c(195,225,255),-0.075,labels=c("July","August","September"),cex=1.5)

legend(x="topleft", legend=c("bud","flower","fruit","seed"), lty=1,lwd=2,
       col=c("yellowgreen","gold1","green4","darkgoldenrod"))
#dev.off()

#2. relationship between DSS and peak phenophase, across all years
sppplot <- "LUAR" #pick species to plot; check 'species' for possibilities
phenoplot <- "flwr" #pick phenophase; check 'phenocats'

#Extract data
plotdat <- outputpars[outputpars[,1]==sppplot&outputpars[,4]==phenoplot,]

#define parameters
DSS <- plotdat[,5]; peakp <- plotdat[,6]

#define color, shape for plotting symbols (shape = transect, color = year)
plotshp <- rep(21, length=dim(plotdat)[1])
siteloc <- substring(plotdat[,3],1,2)
plotshp[siteloc[]=="GB"] <- 22
plotcol <- rep("royalblue4", length=dim(plotdat)[1])
plotcol[plotdat[,2]==2014] <- "lightblue"
plotcol[plotdat[,2]==2015] <- "orange"
plotcol[plotdat[,2]==2016] <- "yellowgreen"
plotcol[plotdat[,2]==2017] <- "mediumorchid2"
plotcol[plotdat[,2]==2018] <- "pink"

#Analysis
DSStest <- lm(peakp~DSS)

#Make a plot: write it as a tiff file or save as pdf
X11(width=6,height=6)
#png(filename="LUARflwr.png",width=6,height=6, units="in", res=1200)
par(omi=c(0,0,0,0), mai=c(0.6,0.6,0.6,0.6), mgp=c(1.45,0.4,0), tck=-0.01,xpd="NA")
plot(DSS, peakp, xlab="Date of Snowmelt", ylab="Peak Flowering",
     pch=plotshp, bg=plotcol, cex=1.5, cex.lab=1.25, xlim=c(90,210),ylim=c(150,240), 
     xaxp=c(90,210,4), yaxp=c(150,240,3), xaxt="n",yaxt="n")

text(c(105,135,165,195),142.5,labels=c("April","May","June","July"), cex=1.25)
text(80, c(165,195,225),labels=c("June","July","August"),srt=90, cex=1.25)


#abline(0,1)
xx <- c(90,210)
yy<-c(coef(DSStest)[1]+coef(DSStest)[2]*90,coef(DSStest)[1]+coef(DSStest)[2]*210)
lines(xx,yy, lty=2)
xx2 <- c(150,210)
yy2 <- c(150,210)
lines(xx2,yy2)
title("Subalpine Lupine Flowering")

#Add a legend
legend(x="topleft", legend=c("RL","GB","2013","2014","2015","2016","2017","2018"),
       pch=c(21,22,21,21,21,21), pt.bg=c("grey","grey","royalblue4","lightblue","orange","yellowgreen","mediumorchid2","pink"))
#dev.off()


#2. relationship between DSS and duration, across all years
sppplot <- "LUAR" #pick species to plot; check 'species' for possibilities
phenoplot <- "flwr" #pick phenophase; check 'phenocats'

#Extract data
plotdat <- outputpars[outputpars[,1]==sppplot&outputpars[,4]==phenoplot,]

#define parameters
DSS <- plotdat[,5]; durp <- plotdat[,7]

#define color, shape for plotting symbols (shape = transect, color = year)
plotshp <- rep(21, length=dim(plotdat)[1])
siteloc <- substring(plotdat[,3],1,2)
plotshp[siteloc[]=="GB"] <- 22
plotcol <- rep("royalblue4", length=dim(plotdat)[1])
plotcol[plotdat[,2]==2014] <- "lightblue"
plotcol[plotdat[,2]==2015] <- "orange"
plotcol[plotdat[,2]==2016] <- "yellowgreen"
plotcol[plotdat[,2]==2017] <- "mediumorchid2"

#Analysis
DSStest <- lm(durp~DSS)

#Make a plot: write it as a tiff file or save as pdf
#tiff(file="SnowPheno.tiff",width=5, height=5, units="in", res=600)
png(filename="LUARflwrdur.png",width=6,height=6, units="in", res=1200)
par(omi=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.5), mgp=c(1.25,0.4,0), tck=-0.01,xpd="NA")
plot(DSS, durp, xlab="Date of Snowmelt", ylab="Duration Flowering (days)",
     pch=plotshp, bg=plotcol, cex=1.75, xlim=c(105,210), 
     xaxp=c(120,210,3), xaxt="n")
text(c(105,135,165,195),1.25, labels=c("May","June","July","August"))



#abline(0,1)
xx <- c(105,210)
yy<-c(coef(DSStest)[1]+coef(DSStest)[2]*105,coef(DSStest)[1]+coef(DSStest)[2]*210)
lines(xx,yy, lty=2)
title("Subalpine Lupine Flowering")

#Add a legend
legend(x="topright", legend=c("RL","GB","2013","2014","2015","2016","2017"),
       pch=c(21,22,21,21,21,21), pt.bg=c("grey","grey","royalblue4","lightblue","orange","yellowgreen","mediumorchid2"))
dev.off()



#3. length of time from flower to seed; or flower to fruit
  #VASI: Bud - fruit; MIAL bud- seed; PEBR bud- fruit; LIGR bud-fruit; 
  #CAMI bud-fruit; ASLE flwr - fruit; ARLA bud-seed; ANOC flwr-seed
  #PLOT: LUAR bud-seed*; POBI - bud - seed*; ERPE bud-seed; CAPA bud-seed?
#4. Peak phenophases across years (especially 2015 vs. other years)
#5. QA / QC for phenophases
#6. Timing of all phenophases for one species (time to reproduction)
#7. Peak wildflower season (in an average year)

###################################
## Additional Analyses & Figures ##
###################################

# maturation affected by SDD & species?
# use first flower to first seed; analysis by species and SDD (per plot)

#First define parameters
species <- unique(PhenoSite[,11]) #first identify unique species
species <- species[order(species)]
years <- unique(PhenoSite[,1])
sites <- unique(PhenoSite[,2])

#Output table
output <- c()

##use for loop to pull out appropriate data for species / plot combos
for(i in 1:length(species)){ #loop for each species
  
  #extract data for that species
  speciesdat<-PhenoSite[PhenoSite[,11]==species[i],]
  
  #loop through each year
  for(j in 1:length(years)){ # loop for each year
    speciesyeardat <- speciesdat[speciesdat[,1]==years[j],]
    
    #break to next if no data for that species / year
    if(dim(speciesyeardat)[1]==0){next}
    
    for(k in 1:length(sites)){ #loop for each phenological site
      #extract data from species dat
      speciesyearsitedat <- speciesyeardat[speciesyeardat[,2]==sites[k],]
      speciesyearsitedat <- speciesyearsitedat[is.na(speciesyearsitedat[,14])==FALSE
                                               & is.na(speciesyearsitedat[,16])==FALSE,]
      
      #stop loop if no data, or if all data is 0 for flowers and or seeds
      if(dim(speciesyearsitedat)[1]==0){next}
      if(sum(speciesyearsitedat[,14])==0|sum(speciesyearsitedat[,16])==0){next}
      
      #first day & last day of observation: find
      days <- speciesyearsitedat[,7]
      day1 <- min(days); daylast <- max(days)     
      
      #first day of flowering: find
      day1f <- min(days[speciesyearsitedat[,14]==1])
      day1s <- max(days[speciesyearsitedat[,16]==1])
      
      #how many observations within 4 days of 1st days? Use as weighting function
      obs1feffort <- length(days[days[]<day1f+4 & days[] > day1f-4])
      obs1seffort <- length(days[days[]<day1s+4 & days[] > day1s-4])
      
      #Pull out SDD
      SDD <- speciesyearsitedat[1,5]
      
      #add to output tables: parameters
      tmpoutput <- c(as.character(species[i]), years[j], as.character(sites[k]), 
                     SDD, day1, daylast, day1f, day1s, obs1feffort, obs1seffort)
      output <- rbind(output, tmpoutput)
    }
  }
}

#Output tables - add dimension names, make appropriate columns numeric
dimnames(output) <- list(c(), c("Species","Year","Site","SDD","firstobs","lastobs", "firstflwr","firstseed","obsflwr","obsseed"))
outputobs <- data.frame(output)
for(i in 4:10){outputobs[,i] <- as.numeric(as.character(outputobs[,i]))}

##Run analysis
daysmat <- outputobs[,8] - outputobs[,7]
spp <- as.factor(outputobs[,1])
year <- as.factor(outputobs[,2])
SDD <- outputobs[,4]
whts <- rowSums(outputobs[,9:10])

test <- lm(daysmat~spp*SDD - 1, weights=whts)



##########################
#Analyze and plot duration of each phenophase vs. DSS & year (four graphs per species)

for(i in 1:length(species)){
  sppdat <- outputpars[outputpars[,1]==species[i],]
  
  #One plot per species
  X11(width=6,height=9)
  par(mfrow=c(4,2), omi=c(0,0,0,0), mai=c(0.4,0.4,0.4,0.4),tck=-0.02, mgp=c(1.25,0.5,0))
  
  
  for(j in 1:length(phenocats)){ #now by phenophase
    #extract data for phenophase
    sppphenodat <- sppdat[sppdat[,4]==names(phenocats)[j],]
    
    #go to next if too few observations; less than 4
    if(dim(sppphenodat)[1]<4){next}
    
    #define parameters
    year <- as.factor(sppphenodat[,2]); DSS <- sppphenodat[,5]; duration <- sppphenodat[,7]
    
    #define color, shape for plotting symbols (shape = transect, color = year)
    plotshp <- rep(21, length=dim(sppphenodat)[1])
    siteloc <- substring(sppphenodat[,3],1,2)
    plotshp[siteloc[]=="GB"] <- 24
    plotcol <- rep("yellowgreen", length=dim(sppphenodat)[1])
    plotcol[sppphenodat[,2]==2014] <- "magenta"
    plotcol[sppphenodat[,2]==2015] <- "orange"
    plotcol[sppphenodat[,2]==2016] <- "navyblue"
    plotcol[sppphenodat[,2]==2017] <- "yellow"
    plotcol[sppphenodat[,2]==2018] <- "purple"
    
    #Analysis 1
    DSStest <- lm(duration~DSS)
    
    #Make a plot
    plot(DSS, duration, xlab="Snowmelt (Julian Days)", ylab=paste("duration(",names(phenocats)[j],")",sep=""),
         pch=plotshp, bg=plotcol, cex=1.5)
    abline(0,1)
    abline(coef(DSStest)[1], coef(DSStest)[2], lty=2)
    title(paste(species[i], names(phenocats)[j], sep="-"))
    
    #Analysis 2
    if(length(unique(year))==1){next}
    yeartest <- lm(duration~year)
    
    #Make a plot
    plot(duration~year, xlab="Year", ylab=paste("duration(",names(phenocats)[j],")",sep=""),
         col=c("yellowgreen","magenta","orange","blue","yellow","purple"))
  }
}


##Check whether length of time (from flowering to fruitset) varies per year, DSS 
#For loop: per species (2 plots)
for(i in 1:length(species)){
  sppdat <- outputpars[outputpars[,1]==species[i],]
  
  #Extract length of time between flowering & seedset
  flwrdat <- sppdat[sppdat[,4]=="flwr",]
  seeddat <- sppdat[sppdat[,4]=="seed",]
  flwrseeddat <- merge(flwrdat,seeddat,by=c("Year","Site"))
  flwrseeddat <- flwrseeddat[,c(1:3,5,6,12)]
  timemat <- flwrseeddat[,6]-flwrseeddat[,5]
  flwrseeddat <- cbind(flwrseeddat,timemat)
  dimnames(flwrseeddat)[[2]] <- c("Year","Site","Species","SDD","peakflower","peakseed","length")
  
  #Define parameters
  year <- as.factor(flwrseeddat[,1])
  SDD <- flwrseeddat[,4]
  maturation <- flwrseeddat[,7]
  if(dim(flwrseeddat)[1]<10){next}
  
  #One plot per species
  X11(width=8,height=5)
  par(mfrow=c(1,2), omi=c(0,0,0,0), mai=c(0.5,0.5,0.5,0.5),tck=-0.02, mgp=c(1.25,0.5,0))
  
  #Analysis; as a function of year
  print(species[i]); print("")
  yeartest <- lm(maturation~year)
  print(summary(yeartest))
  #print(anova(yeartest))
  plot(maturation~year, col=c("yellowgreen","magenta","orange","blue","yellow","purple"), ylab="Days to Seed Dispersal")
  title(species[i])
  
  #Analysis, as a function of SDD
  SDDtest <- lm(maturation ~ SDD)
  print(summary(SDDtest))
  
  #define color, shape for plotting symbols (shape = transect, color = year)
  plotshp <- rep(21, length=dim(flwrseeddat)[1])
  siteloc <- substring(flwrseeddat[,2],1,2)
  plotshp[siteloc[]=="GB"] <- 22
  plotcol <- rep("yellowgreen", length=dim(flwrseeddat)[1])
  plotcol[flwrseeddat[,1]==2014] <- "magenta"
  plotcol[flwrseeddat[,1]==2015] <- "orange"
  plotcol[flwrseeddat[,1]==2016] <- "navyblue"
  plotcol[flwrseeddat[,1]==2017] <- "yellow"
  plotcol[flwrseeddat[,1]==2018] <- "purple"
  
  #Now plot
  plot(SDD, maturation, pch=plotshp, bg=plotcol, xlab="Snowmelt (Julian Days)", ylab="Days to Seed", cex=1.5)
}

