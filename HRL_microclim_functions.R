##Functions for formatting and processing microclimate data.
##Author: Ian Breckheimer
##Created: 15 December 2016

####Function to format an individual HOBO or ibutton file and extract metadata.####
format_micro_csv <- function(csv_name) {
  
  ## Checks to make sure that the file exists in the current working directory.
  stopifnot(file.exists(csv_name))
  print(paste("Now processing file ",csv_name))
  
  ## Reads in the beginning of the file to check the formatting.
  header <- scan(csv_name,nlines=30,what='raw',quiet=TRUE) 
  
  is.formatted <- tolower(header[1]) %in% c("year,month,day,hour,temp,light",
                                            "year,month,day,hour,temperature",
                                            "year,month,day,hour,temp",
                                            "year,month,day,hour,temp,",
                                            "year,month,day,hour,temp,,",
                                            "year,month,day,hour,temp,,,",
                                            "year,month,day,hour,temp,,,,",
                                            "year,month,day ,hour,temp",
                                            "year,month,day,hour,value",
                                            "year,month,day,hour,temp,light,,",
                                            "year,month,day,hour,temp,light,,,",
                                            "year,month,day,hour,temp,light,,,,",
                                            "year,month,day,hour,temp,light,,,,,,,",
                                            "year,month,day,hour,temp,intensity",
                                            "year,month,day,hour,temp,temp")
  
  is.datetime <- unlist(strsplit(header[1],split=","))[1] == "Date/Time"
  is.ibutton <- any(grepl("iButton", header))  # determine whether the sensor is an ibutton or not
  is.hobo <- any(grepl("Plot", header)) && any(grepl("Time", header)) && any(grepl("Title", header))
  if (is.formatted == TRUE){
    file <- read.table(csv_name,skip=1,sep=",")
    DateTime <- paste(file[,2],"/",file[,3],"/",file[,1]," ",file[,4],":00",sep="")
    file <- cbind(DateTime,file[-c(1:4)],NA)
  }
  else if (is.datetime) {
    table <- read.table(csv_name,skip=1,header=FALSE,sep=",")  # remove rows from the top of the data frame
    DateTime <- strptime(table[,1],format="%m/%d/%y %H:%M")
    temp <- as.vector(table[,2])
    file <- data.frame(DateTime,temp,light=NA)
  }
  else if (is.ibutton == TRUE) {   
    # do this if the sensor is an ibutton do this if the ibutton is model DS1921G
    if (any(grepl("DS1921G", header)) == TRUE) {
      file <- read.table(csv_name,skip=15,header=FALSE,sep=",")  # remove rows from the top of the data frame
    }
    if (any(grepl("DS1922L",header)) == TRUE) {
      # do this if the ibutton is model DS1922L
      file <- read.table(csv_name,skip=22,header=FALSE,sep=",")  # remove rows from the top of the data frame
    }
    if (any(grepl("DS2422",header)) == TRUE) {
      # do this if the ibutton is model DS1922
      file <- read.table(csv_name,skip=20,header=FALSE,sep=",")  # remove rows from the top of the data frame
    }
    if (any(grepl("DS1923",header)) == TRUE) {
      # do this if the ibutton is model DS1923
      file <- read.table(csv_name,skip=20,header=FALSE,sep=",")  # remove rows from the top of the data frame
    }
    # do this if there is no header.
    if (class(file) == "factor") {
      if (file[1] == "Date/Time") {
        file = file[-c(1:3)]
      }
      
      if (file[1] == "Unit") {
        file = file[-c(1:2)]
      }
      
      if (file[1] == "Value") {
        file = file[-1]
      }
      
      file = as.character(file)
      tmp.file <- matrix(data = file, ncol = 3, nrow = length(file)/3, byrow = TRUE)  # take data (which is currently in vector format) and convert to matrix format
      file <- as.data.frame(tmp.file)  # convert to data frame
      file <- file[, -2]  # remove the second column (unit)
    }
    if (dim(file)[2] == 1) {
      file = as.character(file)
      tmp.file <- matrix(data = file, ncol = 3, nrow = length(file)/3, byrow = TRUE)  # take data (which is currently in vector format) and convert to matrix format
      file <- as.data.frame(tmp.file)  # convert to data frame
      file <- file[, -2]  # remove the second column (unit)
    }
    if (dim(file)[2] == 2) {
      file = file
    }
    if (dim(file)[2] >= 3) {
      file <- file[, -2]  # remove the second column (unit), and any other columns
    }
    
    light.fill <- array("", dim = c(dim(file)[1], 1))  # make a vector of NAs to fill the light column (ibuttons don't record light)
    file[, 3] = light.fill  # add the vector to the data frame
    
  }  # end ibutton-specific procedure
  else if (is.hobo==TRUE) {
    # do this if the sensor is NOT an ibutton, and is a HOBO
    sample <- header[c((length(header)-10):length(header))]
    sample_split <- unlist(strsplit(sample,split=","))
    is.comma.del <- length(sample)!=length(sample_split)
    if(is.comma.del) {
      file <- read.table(csv_name, header = F, sep = ",",skip=2,fill=TRUE)  # read in the raw data .csv file (as comma delimited) as a dataframe
    }
    else{
      file <- read.table(csv_name,header = F, sep="\t",skip=2,fill=TRUE)    # read in the raw data .csv file (as tab delimited) as a dataframe
    }
    file <- file[,2:3]  # remove all columns except 2 and 3 (date/time and temperature)
  }  # end HOBO-specific procedure
  else{
    print(paste("Could not recognize the formatting of ",csv_name))
  }
  
  names(file) <- c("DateTime", "Temperature")  # name the columns
  
  # Attempts to extract the time-zone from the header.
  if(any(grepl("GMT-07:00", header,fixed=T)) | any(grepl("PDT", header,fixed=T))){
    tz <- "PDT"
  }
  else if(any(grepl("GMT-00:00", header,fixed=T))){
    tz <- "GMT"
  }
  else{
    tz <- "UNK"
  }
  
  # Tests for a valid 4-digit year
  valid_years <- 2007:2050
  datestring <- strsplit(as.character(file$DateTime[1]),split=" ")[[1]][1]
  yeartest <- function(x){any(grepl(as.character(x),datestring))}
  valid.year <- any(sapply(valid_years,FUN=yeartest))
  am.pm <- any(grepl("PM",file$DateTime[1:20]))
  
  # Converts date vector to separate columns for year,month,day,and hour.
  if(class(file$DateTime)[1]=="POSIXct"){
    dateTime <- file$DateTime
  }else if(am.pm && valid.year) {
    dateTime <- strptime(file$DateTime, "%m/%d/%Y %r")
  }else if(am.pm==FALSE && valid.year==TRUE) {
    dateTime <- strptime(file$DateTime, "%m/%d/%Y %H:%M")
  }else if(am.pm==TRUE && valid.year==FALSE) {
    dateTime <- strptime(file$DateTime, "%m/%d/%y %r")
  }else{
    dateTime <- strptime(file$DateTime, "%m/%d/%y %H:%M")
  }
  
  YEAR <- as.numeric(strftime(dateTime,format="%Y"))
  MONTH <- as.numeric(strftime(dateTime,format="%m"))
  DAY <- as.numeric(strftime(dateTime,format="%d"))
  HOUR <- as.numeric(strftime(dateTime,format="%H"))
  MIN <- as.numeric(strftime(dateTime,format="%M"))
  TEMP <- as.numeric(file$Temperature)
  file <- cbind(YEAR,MONTH,DAY,HOUR,MIN,TEMP)
  
  # Extracts the range of dates in the file.
  date_min <- min(dateTime,na.rm=T)
  date_max <- max(dateTime,na.rm=T)
  date_lab <- paste(strftime(date_min,format="%Y-%m"),strftime(date_max,format="%Y-%m"),sep="-")
  
  # Measures the logging interval
  log_int <- dateTime[2] - dateTime[1]
  
  # Extracts the minimum and maximum temperature.
  temp_max <- max(TEMP,na.rm=T)
  temp_min <- min(TEMP,na.rm=T)
  
  # Extracts the number of measurements.
  n_measurements <- length(na.omit(TEMP))
  
  file_attrib <- list(data=file,
                      filename=strsplit(csv_name,split=".csv")[[1]],
                      filepath=paste(getwd(),csv_name,sep="/"),
                      n_measurements=n_measurements,
                      log_interval=log_int,
                      date_min=date_min,
                      date_max=date_max,
                      date_lab=date_lab,
                      temp_min=temp_min,
                      temp_max=temp_max,
                      hobo=is.hobo,
                      ibutton=is.ibutton,
                      formatted=is.formatted,
                      datetime=is.datetime,
                      tz=tz)
  
  return(file_attrib)
}

####Function to batch-format mixed ibutton/HOBO microclimate files in a bunch of directories.####
batch_format_micro_csv <- function(input_paths=getwd(),output_path,file_prefixes,
                                   output_metadata_filename="metadata.txt",overwrite=FALSE){
  
  ##Checks inputs
  stopifnot(length(input_paths)==length(file_prefixes))
  stopifnot(length(output_metadata_filename)==1)
  stopifnot(all(dir.exists(input_paths)))
  stopifnot(dir.exists(output_path))
  stopifnot(is.logical(overwrite) & length(overwrite)==1)
  
  ##Deletes the existing metadata file if it exists
  setwd(output_path)
  if(file.exists(output_metadata_filename)){
    file.remove(output_metadata_filename)
  }
  
  ##Checks to see how many input files there are.
  csv_all <- c()
  for (i in 1:length(input_paths)){
    csv_files <- list.files(input_paths[i],pattern=".csv$")
    csv_all <- c(csv_all,csv_files)
  }
  nfiles <- length(csv_all)
  print(paste("Now processing ",nfiles," microclimate files."))
  
  ##Sets up a file counter.
  file_n <- 0
  
  for (i in 1:length(input_paths)){
    setwd(input_paths[i])
    flush.console()
    print(paste("Now processing files in folder ",input_paths[i]))
    files <- list.files(".",pattern=".csv$") # lists all the .csv files in the working directory (!!! WARNING: DON'T INCLUDE '.csv' IN ANY FILE NAMES IN THE WORKING DIRECTORY, '.csv' SHOULD ONLY APPEAR AS A FILE EXTENSION !!!). This script thinks that any file with the string '.csv' in the file name is a .csv file
    dfs <- lapply(X = files, FUN = format_micro_csv)  # process all .csv files in the working directory into a list of lists.
    
    ## write the modified files as .csv files to an output folder
    
    # set the working directory to the folder that will collect the output files
    setwd(output_path)
    outfiles <- length(list.files())
    if(outfiles>0){
      print(paste("Output folder already contains ", outfiles,"files."))
    }
    
    for (j in 1:length(dfs)) {
      datavals <- dfs[[j]]$data
      in_name <- dfs[[j]]$filename
      out_name <- paste(paste(file_prefixes[i],in_name,dfs[[j]]$date_lab,sep="_"),".csv",sep="")
      meta <- data.frame(out_filename=out_name,dfs[[j]][-1])
      if(overwrite==TRUE & file.exists(out_name)){
        print(paste("Overwriting file", out_name))
        write.csv(datavals, file = out_name, row.names = FALSE)
      }else if(overwrite==FALSE & file.exists(out_name)){
        print(paste("Skipping existing file", out_name))
      }else{
        print(paste("Writing file",out_name))
        write.csv(datavals, file = out_name, row.names = FALSE)
      }
      
      # Writes the metadata, appending rows to an existing file if it already exists.
      if (length(list.files(".",pattern=output_metadata_filename)) == 0){
        write.table(meta,file=output_metadata_filename, sep=",",row.names = FALSE, append = FALSE) 
      }
      else{
        write.table(meta,file=output_metadata_filename, sep=",",row.names = FALSE, 
                    col.names = FALSE, append = TRUE)
      }
      file_n <- file_n+1
      print(paste("Completed processing file ",file_n," of ",nfiles))
    }
    
  }
  
  ##Returns output metadata
  setwd(output_path)
  meta_output <- read.csv(output_metadata_filename)
  return(meta_output)
}

####Function to estimate snow cover duration from a series of .csv microclimate files.####
batch_extract_snow_vars <- function(input_path,input_metadata_filename="metadata.txt",
                                    figure_path,output_path,output_metadata_filename,
                                    range_threshold=1,max_threshold=2,overwrite=FALSE){
  ##Loads required packages
  require(data.table)
  
  ##Checks inputs
  stopifnot(dir.exists(input_path))
  stopifnot(dir.exists(figure_path))
  stopifnot(dir.exists(output_path))
  
  setwd(input_path)
  stopifnot(file.exists(input_metadata_filename))
  stopifnot(!file.exists(output_metadata_filename)|overwrite==TRUE)
  
  files <- list.files(".",pattern=".csv$")  # the data files for each temperature sensor to be analyzed
  meta <- read.table(input_metadata_filename,sep=",",header=TRUE)
  
  calibration <- c()
  calibration.type <- c()
  stand <- c()
  plot <- c()
  year <- c()
  snow_appearance_date <- c()
  snow_disappearance_date <- c()
  snow_cover_duration <- c()  # in days
  minimum_soil_temp <- c()
  
  # start the clock to record how long the code takes to run
  start_t <- Sys.time()
  
  nfiles <- length(files)
  for (k in 1:nfiles) {
    
    ##Makes sure we are in the right directory
    setwd(input_path)
    
    ##Prints progress.
    flush.console()
    print(paste("Now Processing file: ",files[k],"(",k," of",nfiles,")"))
    
    ############################################## READING IN DATA FROM ONE FILE
    
    d <- read.csv(files[k])
    d <- d[complete.cases(d[, 1:5]), ]
    
    # Name the columns
    names(d)[1] <- "YEAR"
    names(d)[2] <- "MONTH"
    names(d)[3] <- "DAY"
    names(d)[4] <- "HOUR"
    names(d)[5] <- "MIN"
    names(d)[6] <- "TEMP"
    
    d$Date <- as.Date(paste(d$MONTH, "/", d$DAY, "/", d$YEAR, sep = ""), format = "%m/%d/%Y")
    d$DOY <- as.numeric(format(d$Date, format = "%j"))  # find the unique days
    d <- data.table(d,key="Date")
    
    ######################################## 
    
    # EXTRACT STAND INFO FROM FILENAME
    
    stand[k] <- strsplit(files[k], "_")[[1]][2]
    plot[k] <- strsplit(files[k], "_")[[1]][3]
    year[k] <- max(d$YEAR)
    
    ############################################### 
    
    # EVALUATE SNOW COVER CRITERIA
    
    # Create an empty data table to store values
    days <- unique(d$Date)
    ndays <- length(days)
    
    daily <- data.table(date=unique(d$Date),
                        range=rep(NA,ndays),
                        mean=rep(NA,ndays),
                        rangethresh=rep(NA,ndays),
                        maxthresh=rep(NA,ndays),
                        snow=rep(NA,ndays))
    setkey(daily,"date")
    
    # Calculate the mean daily temperature and temperature range for each day:
    daily[[2]] <- d[,diff(range(TEMP)),by=Date][[2]]
    daily[[3]] <- d[,mean(TEMP,na.rm=T),by=Date][[2]]
    daily[[4]] <- d[,diff(range(TEMP)) < range_threshold,by=Date][[2]]
    daily[[5]] <- d[,max(TEMP) < max_threshold,by=Date][[2]]
    
    ############################################## 
    
    # DETERMINE CALIBRATION TEMPERATURE
    calibration.temp <- mean(daily$mean[daily$rangethresh == TRUE])  # calculate the mean temp for all the days that the temp didn't exceed range_threshold
    calibration[k] <- ifelse(!is.na(calibration.temp), yes = calibration.temp, no = 0)  # If NA, set calibration = 0
    d$TEMP.calib <- d$TEMP - calibration[k]  # recalibrate data 
    
    ################################################ 
    
    ## EVALUATE WHETHER OR NOT SNOW COVERED THE SENSOR ON EACH DATE BASED ON THE ABOVE TWO CRITERIA
    
    daily[[6]] <- daily[[4]] & daily[[5]]  # will store the algorithm's evaluation of snow cover for each date (1=snow present, 0=snow absent)
    d.snow <- subset(daily,snow==TRUE)
    ################################################ 
    
    ## SUMMARIZING
    snow_appearance_date[k] <- as.character(min(d.snow$date))  # first day when snow covered sensor
    snow_disappearance_date[k] <- as.character(max(d.snow$date))  # last day when snow covered sensor
    snow_cover_duration[k] <- sum(d.snow$snow)  #'snow cover duration' the total number of days with snow cover
    minimum_soil_temp[k] <- min(d$TEMP,na.rm=TRUE) #Winter minimum soil temperature
    
    
    ## PLOT SOIL TEMPERATURE AND THE SNOW COVER ALGORITHM OUTPUT TO MAKE SURE OUTPUT IS REASONABLE
    
    # Save figure as pdf
    setwd(figure_path)
    population <- strsplit(files[k], ".csv")[[1]]
    graph.file <- paste(population, ".pdf", sep = "")
    pdf(file = graph.file, width = 10, height = 7)
    
    # Left Y axis
    par(mar = c(4, 6, 3, 6))
    plot(d$Date, d$TEMP.calib, main = population, axes = F, xlab = "", ylab = "", pch = 20, col = "red")
    axis(2, col = "red", col.axis = "red", col.ticks = "red")
    leftY <- expression(paste("Temperature", degree, "C"))
    text(par("usr")[1] - 30, par("usr")[3] + ((par("usr")[3] + par("usr")[4])/2), adj = 0.5, leftY, srt = 90, xpd = TRUE, 
         col = "red")
    
    # Right Y axis
    par(new = TRUE)
    plot(daily$date, daily$snow, pch = 20, col = "blue", ylim = c(0, 1.2), axes = FALSE, xaxt = "n", yaxt = "n", 
         xlab = "", ylab = "")
    axis(4, col = "blue", col.axis = "blue", col.ticks = "blue", yaxp = c(0, 1, 1), labels = F)
    rightY <- paste("Snow cover")
    text(par("usr")[2] + 30, par("usr")[3] + ((par("usr")[3] + par("usr")[4])/2), adj = 0.5, rightY, srt = 270, xpd = TRUE, 
         col = "blue")
    text(par("usr")[2] + 15, 0, 0, srt = 270, xpd = T, col = "blue")
    text(par("usr")[2] + 15, 1, 1, srt = 270, xpd = T, col = "blue")
    
    # X axis
    axis.Date(1, d$Date, at = seq(from = min(d$Date), to = max(d$Date), by = "months"))
    
    dev.off()
    
    ################################################ 
    
  }
  
  ##Computes elapsed time.
  stop_t <- Sys.time()
  print(paste("Run took ",stop_t - start_t))
  
  # Consolidate summarized results for each sensor into one data frame
  output <- data.frame(files, 
                       year, 
                       stand, 
                       plot, 
                       calibration, 
                       snow_appearance_date, 
                       snow_disappearance_date,
                       snow_cover_duration,
                       minimum_soil_temp)
  
  ##Match the snow cover information to the sensor metadata.
  out_merged <- unique(merge(meta,output,by.x="out_filename",by.y="files",all.x=T))
  
  ##Data quality flags:
  end_doy <- as.numeric(format(as.Date(out_merged$date_max),format="%j"))
  out_merged$flag_sensor_fail <- ((as.Date(out_merged$snow_disappearance_date) - as.Date(out_merged$date_max)) >= -1) & end_doy < 152 
  out_merged$flag_temp_high <- out_merged$temp_max > 90
  out_merged$flag_temp_low <- out_merged$temp_min < -30
  out_merged$flag_high_calib <- out_merged$calibration >= 2 | out_merged$calibration <= -2
  out_merged$flag_no_snow <- out_merged$snow_cover_duration <= 14
  out_merged$flag_short_record <- as.Date(out_merged$date_max) -  as.Date(out_merged$date_min) < 100
  out_merged$flagged <- with(out_merged, flag_sensor_fail | flag_temp_high | flag_temp_low | flag_high_calib | flag_no_snow | flag_short_record)
  
  ##Moves .csv and .pdf graphics of flagged files to new directories
  out_flagged <- out_merged[out_merged$flagged==TRUE,]
  flagged_csvs <- out_flagged$out_filename
  setwd(output_path)
  dir.create("./flagged")
  file.copy(paste(input_path,"/",as.character(flagged_csvs),sep=""),"./flagged")
  flagged_pdfs <- sub(".csv",".pdf",out_flagged$out_filename)
  setwd(figure_path)
  dir.create("./flagged")
  file.copy(flagged_pdfs,"./flagged")
  
  ##Moves .csv and .pdf graphics of unflagged files to a new directory
  out_unflagged <- out_merged[out_merged$flagged==FALSE,]
  unflagged_csvs <- out_unflagged$out_filename
  setwd(output_path)
  dir.create("./unflagged")
  file.copy(paste(input_path,"/",as.character(unflagged_csvs),sep=""),"./unflagged")
  unflagged_pdfs <- sub(".csv",".pdf",out_unflagged$out_filename)
  setwd(figure_path)
  dir.create("./unflagged")
  file.copy(unflagged_pdfs,"./unflagged")
  
  ##subsets data for the unflagged files
  out_unflagged<- out_merged[out_merged$flagged==FALSE,]
  
  # Save output file
  setwd(output_path)
  write.table(out_unflagged, file = output_metadata_filename, sep = ",", row.names = FALSE)
  return(out_merged)
}

batch_clean_air_temps <- function(input_path,input_metadata_filename,
                                      output_path,output_metadata_filename,
                                      figure_path,
                                      guess_tz="Etc/GMT-7",
                                      temp_spike_thresh=10,
                                      min_temp_thresh=-20,
                                      max_temp_thresh=50,
                                      max_temp_hr=17,
                                      overwrite=TRUE){
  ##Sets up workspace
  require(xts)
  require(psych)
  Sys.setenv(TZ=guess_tz)
  
  ##Checks inputs
  stopifnot(length(input_metadata_filename)==1)
  stopifnot(length(output_metadata_filename)==1)
  stopifnot(dir.exists(input_path))
  stopifnot(dir.exists(output_path))
  stopifnot(dir.exists(figure_path))
  stopifnot(is.logical(overwrite) & length(overwrite)==1)
  setwd(input_path)
  stopifnot(file.exists(input_metadata_filename))
  
  ##Checks to see how many input files there are.
  airtemp_files <- list.files(input_path,pattern=".csv$")
  nfiles <- length(airtemp_files)
  
  ##Creates folders for output
  setwd(output_path)
  if(!dir.exists("./flagged")){
    dir.create("./flagged")
  }
  if(!dir.exists("./unflagged")){
    dir.create("./unflagged")
  }
  
  ##Checks to make sure we've got the same number files as we do metadata lines
  setwd(input_path)
  air_meta <- read.table(input_metadata_filename,sep=",",header=TRUE)
  stopifnot(nrow(air_meta)==nfiles)
  
  print(paste("Now processing ",nfiles," air-temp files."))
  
  ##Sets up a file counter.
  file_n <- 0
  
  ##Creates empty data frame for metadata.
  nfiles <- length(airtemp_files)
  metadata <- data.frame(filename=rep("NA",nfiles),
                         cleaned_min=rep(NA,nfiles),
                         cleaned_max=rep(NA,nfiles),
                         mean_hr_max=rep(NA,nfiles),
                         flag_missing_data=rep(NA,nfiles),
                         flag_wrong_tz=rep(NA,nfiles))
  
  metadata$filename <- as.character(metadata$filename)
  
  for (i in 1:length(airtemp_files)){
    
    # Ensures we are looking at the right directory
    setwd(input_path)
    
    # Prints progress to console.
    flush.console()
    print(paste("Now processing ",airtemp_files[i],". File ",i," of ",
                length(airtemp_files)))
    #Reads data and converts to date and time-series format.
    data <- read.csv(airtemp_files[i],header=TRUE)
    data <- data[complete.cases(data),]
    data$datestring <- paste(data$YEAR,data$MONTH,data$DAY,sep="-")
    data$timestring <- paste(data$HOUR,data$MIN,"00",sep=":")
    data$datetime <- strptime(paste(data$datestring,data$timestring),
                              format="%Y-%m-%d %H:%M:%S",tz=guess_tz)
    tempts <- xts(data$TEMP,order.by=data$datetime)
    tzone(tempts) <- "Etc/GMT-7"
    
    # Save figure as pdf
    setwd(figure_path)
    plotname <- strsplit(airtemp_files[i], ".csv")[[1]]
    figpath <- paste(figure_path,plotname,".pdf", sep = "")
    pdf(file = figpath, width = 10, height = 7)
    
    # Left Y axis
    par(mar = c(4, 6, 3, 6))
    plot(data$datetime,data$TEMP,type='l', main = plotname, xlab='Date', 
         ylab=expression(paste("Temperature", degree, "C")),col='blue',
         axes = T,ylim=c(-10,35),pch = 20)  
    dev.off()
    
    ## Checks for and eliminates values outside of the physical range 
    ## and temperature spikes that could be caused by direct sun.
    tempts[tempts < min_temp_thresh] <- NA
    tempts[tempts > max_temp_thresh] <- NA
    tlagged <- lag.xts(tempts,k=1)
    delta <- tempts-tlagged
    tempts[delta > temp_spike_thresh] <- NA
    
    ##Calculates the circadian mean hour of the day with the maximum temperature.
    hr_max_ts <- do.call(rbind, lapply(split(tempts,"days"), function(x) x[which.max(x)]))
    hrs <- as.numeric(format(index(hr_max_ts),format="%H"))
    mean_hr_max <- circadian.mean(hrs,hours=TRUE)
    
    ##Flags time-series with problems
    if(sum(is.na(tempts)) > 10){
      flag_missing_values <- TRUE
    }else{
      flag_missing_values <- FALSE
    }
    
    if((mean_hr_max-max_temp_hr) > 4 | (mean_hr_max-max_temp_hr) < -4){
      flag_wrong_tz <- TRUE
    }else{
      flag_wrong_tz <- FALSE
    }
    outts <- data.frame(DATE=index(tempts),TZ=air_meta$tz[i],TEMP=tempts)
    setwd(output_path)
    if(flag_missing_values==FALSE & flag_wrong_tz==FALSE){
      write.csv(outts,paste("./unflagged/",airtemp_files[i]),row.names=FALSE)
    }else{
      write.csv(outts,paste("./flagged/",airtemp_files[i]),row.names=FALSE)
    }
    
    ##Assembles metadata.
    metadata$filename[i] <- airtemp_files[i]
    
    ## Gets the cleaned temp range
    metadata$cleaned_min[i] <- min(tempts,na.rm=TRUE)
    metadata$cleaned_max[i] <- max(tempts,na.rm=TRUE)
    
    ## Gets the hour of mean temperature.
    metadata$mean_hr_max[i] <- mean_hr_max
    
    metadata$flag_missing_data[i] <- flag_missing_values 
    metadata$flag_wrong_tz[i] <- flag_wrong_tz
  }
  
  metadata_unflagged <- metadata[metadata$flag_missing_data==FALSE &
                                 metadata$flag_wrong_tz==FALSE,]
  
  meta_all <- merge(air_meta,metadata_unflagged,by.x="out_filename",by.y="filename",all.x=T)
  
  ##Writes cleaned unflagged files to output.
  setwd(output_path)
  write.table(meta_all,output_metadata_filename,sep="\t",row.names=FALSE)
  return(metadata)
}

