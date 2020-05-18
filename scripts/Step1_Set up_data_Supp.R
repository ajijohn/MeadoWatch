## Author: Meera Lee Sethi
## Last worked on: 12-27-2019
## Script for prepping the dataset used in "Early snowmelt and warmer, drier summers shrink post-flowering transition times in subalpine wildflowers"

#Sets up workspace
library(dplyr)
library(tidyr)
library(reshape)
library(reshape2)
library(runjags)
library(data.table)

#####
## Reads in phenodata and preps it for modeling:
phenodat <- read.csv(file="./Phenology_Env_2010to2015_gapfilled.csv", header=TRUE)
phenodat[,12:23][is.na(phenodat[,12:23])] <- 0 # Converts NAs in observation columns to 0s
phenodat <- mutate(phenodat, Flwr = ifelse(PercFlwr>0, 1, 0)) #Converts percentages to presence-absence data
phenodat <- mutate(phenodat, Fruit = ifelse(PercFruit>0, 1, 0))
phenodat <- mutate(phenodat, Disp = ifelse(PercDispFruit>0, 1, 0))

# Removes species that appear in less than 4 plot-year combinations or where seeds were observed in fewer than 30% of all plot-year combinations
pheno_brief <- filter(pheno_brief, Species!="Epilobium alpinum")
pheno_brief <- filter(pheno_brief, Species!="Caltha leptosepala")
pheno_brief <- filter(pheno_brief, Species!="Cassiope mertensiana")
pheno_brief <- filter(pheno_brief, Species!="Erythronium grandifloru")
pheno_brief <- filter(pheno_brief, Species!="Gentiana calycosa")
pheno_brief <- filter(pheno_brief, Species!="Lomatium martindalei")
pheno_brief <- filter(pheno_brief, Species!="Mitella breweri")
pheno_brief <- filter(pheno_brief, Species!="Pedicularis groendlandi")
pheno_brief <- filter(pheno_brief, Species!="Pedicularis ornithorhyn") 
pheno_brief <- filter(pheno_brief, Species!="Penstemon confertus") 
pheno_brief <- filter(pheno_brief, Species!="Phlox diffusa")
pheno_brief <- filter(pheno_brief, Species!="Phyllodoce glanduliflor") 
pheno_brief <- filter(pheno_brief, Species!="Polemonium californicum")
pheno_brief <- filter(pheno_brief, Species!="Rhododendron albiflorum") 
pheno_brief <- filter(pheno_brief, Species!="Rubus lasiococcus") 
pheno_brief <- filter(pheno_brief, Species!="Salix barclayi")
pheno_brief <- filter(pheno_brief, Species!="Sorbus sitchensis")
pheno_brief <- filter(pheno_brief, Species!="Spiraea splendens")
pheno_brief <- filter(pheno_brief, Species!="Xerophyllum tenax")
pheno_brief <- filter(pheno_brief, Species!="Claytonia sibirica")

pheno_brief <- rename(pheno_brief, c(Species_num="species")) #Standardizes species column name 

pheno_tidy <- pheno_brief %>% #Tidies the data
  gather(phenophase, presence, Flwr:Disp)
pheno_tidy <- pheno_tidy%>% arrange(Year, Site,Species, DOY)

gdd <- read.csv(file="./th_para_sites_gdd_1_2009_10_2015.csv", header=TRUE) ## Reads in gdd data:

gdd <- melt(gdd,id.vars=c("Date","DOY","Year")) #Formats dates in gdd data 
levels(gdd$variable) <- c("1490","1570","1680","1791","1901") # adds elevation identifiers
gdd<- plyr::rename(gdd, c("variable"="elev"))


phenotidy_plots <- unique(data.frame(pheno_tidy$Site, #Creates new dataset with plot info only
                                     pheno_tidy$Year,
                                     pheno_tidy$Elevation,
                                     pheno_tidy$SDDDOY,
                                     pheno_tidy$snow_cover_duration))
colnames(phenotidy_plots) <- c("Site","Year", "Elevation","SDDDOY", "snowcover")

phenotidy_plots$gddsdd <- NA #Adds columns for cumulative gdd25 and gdd 50
phenotidy_plots$gdd25 <- NA
phenotidy_plots$gdd50 <- NA

nrow <- nrow(phenotidy_plots) # Computes cumulative gdd at snowmelt day and 25 days later
for(i in 1:nrow){
  phenotidy_plots$gddsdd[i] <- base::subset(gdd$value, gdd$elev==phenotidy_plots$Elevation[i] & gdd$Year==phenotidy_plots$Year[i] & gdd$DOY==phenotidy_plots$SDDDOY[i])
  phenotidy_plots$gdd25[i] <- base::subset(gdd$value, gdd$elev==phenotidy_plots$Elevation[i] & gdd$Year==phenotidy_plots$Year[i] & gdd$DOY==(phenotidy_plots$SDDDOY[i]+25))
  phenotidy_plots$gdd50[i] <- base::subset(gdd$value, gdd$elev==phenotidy_plots$Elevation[i] & gdd$Year==phenotidy_plots$Year[i] & gdd$DOY==(phenotidy_plots$SDDDOY[i]+50))
}

# computes gdd25 and gdd 50 
phenotidy_plots$gdd50 <- as.numeric(round(phenotidy_plots$gdd50)-as.numeric(phenotidy_plots$gdd25))
phenotidy_plots$gdd25 <- as.numeric(round(phenotidy_plots$gdd25)-as.numeric(phenotidy_plots$gddsdd))


## Merges gdd data back into tidy data
pheno_tidy_nofruit <- merge(pheno_tidy, phenotidy_plots, by =c("Site","Year", "Elevation","SDDDOY"))


### Setting up data to run in the model
setDT(pheno_tidy_nofruit)[, phase := .GRP, by = phenophase] #maps phenophases to unique numerical indices
setDT(pheno_tidy_nofruit)[, plot := .GRP, by = Site] #maps plots to unique numerical indices
setDT(pheno_tidy_nofruit)[, year := .GRP, by = Year] #maps year to unique numerical indices
setDT(pheno_tidy_nofruit)[, species := .GRP, by = Species] #maps species to unique numerical indices
species_ids <- distinct(pheno_tidy_nofruit[, c(8, 23)]) #table of species with id numbers

# Centers and scales predictors
pheno_tidy_nofruit <- droplevels(pheno_tidy_nofruit)
pheno_tidy_nofruit$scaledDOY <- scale(pheno_tidy_nofruit$DOY)
pheno_tidy_nofruit$scaledSDDDOY <- scale(pheno_tidy_nofruit$SDDDOY)
pheno_tidy_nofruit$scaledgdd25 <- scale(pheno_tidy_nofruit$gdd25)
pheno_tidy_nofruit$scaledgdd50 <- scale(pheno_tidy_nofruit$gdd50)
pheno_tidy_nofruit$scaledsoilmoist <- scale(pheno_tidy_nofruit$moist_days)

# Separates data into training and testing subsets

set.seed(42)
pheno_tidy_nofruit<- tibble::rownames_to_column(pheno_tidy_nofruit) #adds rownames as a column
train_data <- pheno_tidy_nofruit%>%
  group_by(Species, Year) %>%
  sample_frac(0.9)
write.csv(train_data,"./final_train_data.csv")
test_data <- anti_join(pheno_tidy_nofruit, train_data,  by="rowname")
write.csv(test_data,"./final_test_data.csv")