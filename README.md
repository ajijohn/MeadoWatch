# MeadoWatch
MeadoWatch Analysis Project

This github repository contains the raw data, cleaned data, scripts (R markdown), figures and output files from analysis characterizing the wildflower season and visitor sentiment. Specifically:

The main folder contains the markdown files (with analyses), including:
1. MW_DataWrangling.RMD: this markdown reads in raw data, and filters out outliers and plot-species combinations with too few observations. This script takes a long time to run, and could probably use some optimizing.
2. MW_SDDModeling.RMD: this markdown reads in SDD data (processed by MW_Microclimate.RMD) from collected from various locations with Mt Rainier NP, and gap fills missing MeadoWatch observations.
3. MW_Modelfitting(byyearspeciesplot).RMD:  this markdown reads in clean phenology data and fits phenological curves to data from each species / year / plot.
4. MW_Modelfitting(byspecies).RMD: this markdown reads in clean phenology data and clean snow disappearance data and fits a phenological model to all data collected from each species (assuming a linear relationship between SDD and peak flowering).
5. MW_PredictedFloweriness.RMD: this markdown reads in species specific parameters (description the relationship between peak flowering and snow for each species), and estimates flowering richness along each trail.
6. MW_WTAanalysis: This markdown reads in WTA data and explores it and predicted flowering (from MeadoWatch data)

The data folder contains
1. Phenology data
2. WTA data
3. Snow disappearance data

The clean data folder contains
1. Phenology data (minus outliers, plot-year-species data with less than a certain number of observations)
2. MeadoWatch data (with predictions for missing data)

The output folder contains
1. parameters from maximum likelihood models (fit per species / year / plot, fit per species)
2. AIC from models
3. outliers

The figs folder contains various figures (many also plotted within the markdown file)
