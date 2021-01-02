##################################
# MeadoWatch WTA Analysis
# Adapted by Aji John #
# original by Alessio Benedetti           #
# server.R file                  #
##################################

library(shiny)
library(tidyverse)
library(leaflet.extras)
library(rvest)

#####################
# SUPPORT FUNCTIONS #
#####################


##################
# DATA WRANGLING #
##################


PhenoSite_clean <- read.csv("cleandata/PhenoSite_Clean.csv", header=TRUE)

MW_PhenoSite_2013_2018 <- read_csv("data/MW_SiteDat_2013_2018.csv")

AllSites_2018 <- MW_PhenoSite_2013_2018 %>% 
  filter(Year == 2018 & Transect %in% c( "Reflection Lakes","Glacier Basin")) %>% 
  as.data.frame()
fp <- read_csv('output/FloweringProbabilities.csv')
meadow_species <- unique(PhenoSite_clean$Species)

fr <- read_csv('output/FloweringRichness.csv')
years <- unique(fp$year)

# tidy & enrich datafram

# support structures

################
# SERVER LOGIC #
################

shinyServer(function(input, output) {
  
  # trails map
  output$trailsMap <- renderLeaflet({
    leaflet(data=AllSites_2018) %>% addProviderTiles(providers$Stamen.Watercolor, group = "Stamen Watercolor", options = providerTileOptions(noWrap = TRUE)) %>%#, minZoom = 4)) %>%
      addProviderTiles(providers$OpenStreetMap.Mapnik, group = "Open Street Map", options = providerTileOptions(noWrap = TRUE)) %>%
      addProviderTiles(providers$NASAGIBS.ViirsEarthAtNight2012, group = "Nasa Earth at Night", options = providerTileOptions(noWrap = TRUE)) %>%
      addProviderTiles(providers$Stamen.TerrainBackground, group = "Stamen Terrain Background", options = providerTileOptions(noWrap = TRUE)) %>%
      addProviderTiles(providers$Esri.WorldImagery, group = "Esri World Imagery", options = providerTileOptions(noWrap = TRUE)) %>%
      addFullscreenControl() %>%
      setView(	-121.7203,46.8523, zoom = 10) %>% 
      addMarkers(
        ~longitude,
        ~latitude,
        icon = makeIcon(
          iconUrl = "flowers_1.png",
          #shadowUrl = "flowers.png",
          shadowAnchorX = -1, shadowAnchorY = -2
        ),
        clusterOptions = markerClusterOptions()
      ) %>%
      addLayersControl(
        baseGroups = c("Stamen Watercolor","Open Street Map","Nasa Earth at Night","Stamen Terrain Background","Esri World Imagery"),
        position = c("topleft"),
        options = layersControlOptions(collapsed = TRUE)
      )
  })
  
  
  # DT table
  output$phenoSiteDataTable <- renderDataTable(
    PhenoSite_clean,
    filter = "top"

  )
  

  # ggplot2 charts
  output$categorySelectComboChart <- renderUI({
    selectInput("selectedCategoryChart","Select a species:", meadow_species)
  })
  
  
  speciesGgplot1 <- reactive(fp[,c('year','trail','DOY',input$selectedCategoryChart)])
  speciesGgplot2 <- reactive(fp[,c('year','trail','DOY',input$selectedCategoryChart)])
  speciesGgplot3 <- reactive(fr)
  speciesGgplot4 <- reactive(fr)
  
  output$ggplot2Group1 <- renderPlot({
    
     
    g1 <- speciesGgplot1() %>% 
      filter(trail=='RL') %>% 
      mutate(yearf = as.factor(year)) %>% 
      ggplot() + geom_point(aes_string(x="DOY",y=input$selectedCategoryChart,color="yearf")) + labs(color='Year' , title="Species - [Refections Lakes Trail]", x ="Day of the year (DOY)", y = paste0("P (flowering ) of ", input$selectedCategoryChart)) +theme_classic() 
    print(g1)
    
  })
  
  output$ggplot2Group2 <- renderPlot({
    
    
    g2 <- speciesGgplot2() %>% 
      filter(trail=='GB') %>% 
      mutate(yearf = as.factor(year)) %>% 
      ggplot() + geom_point(aes_string(x="DOY",y=input$selectedCategoryChart,color="yearf")) + labs(color='Year',title="Species - [Glacier Basin]", x ="Day of the year (DOY)", y = paste0("P (flowering ) of ", input$selectedCategoryChart)) +theme_classic() 
    print(g2)
    
  })
  
  
  # selections for richness
  output$yearSelectCombo <- renderUI({
    selectInput("yearsCombo","Select a year:", years)
  })
  
  output$speciesCombo <- renderUI({
    selectInput("speciesCombo","Select a species:", meadow_species)
  })
  
  
  output$ggplot2Group3 <- renderPlot({
    
    
    g3 <- speciesGgplot3() %>% 
      filter(trail=='RL') %>% 
      mutate(yearf = as.factor(year)) %>% 
      ggplot() + geom_smooth(aes_string(x="DOY",y="rich",color="yearf")) + labs(color='Year' , title="Richness - [Refections Lakes Trail]", x ="Day of the year (DOY)", y = paste0(" ", "")) +theme_classic() 
    print(g3)
    
  })
  
  output$ggplot2Group4 <- renderPlot({
    
    
    g4 <- speciesGgplot4() %>% 
      filter(trail=='GB') %>% 
      mutate(yearf = as.factor(year)) %>% 
      ggplot() + geom_smooth(aes_string(x="DOY",y="rich",color="yearf")) + labs(color='Year',title="Richness - [Glacier Basin]", x ="Day of the year (DOY)", y = paste0("", "")) +theme_classic() 
    print(g4)
    
  })  
  
})