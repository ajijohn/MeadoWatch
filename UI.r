##################################
# MeadoWatch WTA Analysis
# Adapted by Aji John #
# original by Alessio Benedetti           #
# ui.R file                      #
##################################

library(leaflet)
library(shinydashboard)
library(collapsibleTree)
library(shinycssloaders)
library(DT)
library(tigris)

###########
# LOAD UI #
###########

shinyUI(fluidPage(
  
  # load custom stylesheet
  includeCSS("www/style.css"),
  
  # load google analytics script
  #tags$head(includeScript("www/google-analytics-bioNPS.js")),
  
  # remove shiny "red" warning messages on GUI
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  ),
  
  # load page layout
  dashboardPage(
    
    skin = "green",
    
    dashboardHeader(title="MeadoWatch/WTA Analysis", titleWidth = 300),
    
    dashboardSidebar(width = 300,
                     sidebarMenu(
                       HTML(paste0(
                         "<br>",
                         "<a><img style = 'display: block; margin-left: auto; margin-right: auto;' src='mountain.svg' width = '186'></a>",
                         "<br>",
                         "<br>"
                       )),
                       menuItem("Home", tabName = "home", icon = icon("home")),
                       menuItem("Hikes", tabName = "map", icon = icon("thumbtack")),
                       menuItem("Phenology at sites", tabName = "table", icon = icon("table")),
                       menuItem("Model fitting curves", tabName = "charts", icon = icon("stats", lib = "glyphicon")),
                       menuItem("Flowering richness", tabName = "richness", icon = icon("map marked alt")),
                       menuItem("MW WTA analysis", tabName = "wta", icon = icon("random", lib = "glyphicon")),
                       menuItem("Releases", tabName = "releases", icon = icon("tasks")),
                       HTML(paste0(
                         "<br><br><br><br><br><br><br><br><br>",
                         "<table style='margin-left:auto; margin-right:auto;'>",
                         "<tr>",
                          "</tr>",
                         "</table>",
                         "<br>"),
                         HTML(paste0(
                           "<script>",
                           "var today = new Date();",
                           "var yyyy = today.getFullYear();",
                           "</script>",
                           "<div>Icons made by <a href='https://www.flaticon.com/authors/freepik' title='Freepik'>Freepik</a></div>",
                           "<p style = 'text-align: center;'><small>Template credit <a>alessiobenedetti.com</a>-<script>document.write(yyyy);</script></small></p>")
                         ))
                     )
                     
    ), # end dashboardSidebar
    
    dashboardBody(
      
      tabItems(
        
        tabItem(tabName = "home",
                
                # home section
                includeMarkdown("www/home.md")
                
        ),
        
        tabItem(tabName = "map",
                
                # trails map section
                leafletOutput("trailsMap") %>% withSpinner(color = "green")
                
        ),
        
        tabItem(
          # pheno sites observations
          tabName = "table", dataTableOutput("phenoSiteDataTable") %>% withSpinner(color = "green")
          
        ),
        
        tabItem(tabName = "wta", 
                
                # include the WTA Analysis
                includeMarkdown("www/WTA.md")

                
        ),
        
        tabItem(tabName = "charts",
                
                # ggplot2 species charts section
                includeMarkdown("www/charts.md"),
                fluidRow(column(3, uiOutput("categorySelectComboChart"))),
                column(6, plotOutput("ggplot2Group1") %>% withSpinner(color = "green")),
                column(6, plotOutput("ggplot2Group2") %>% withSpinner(color = "green"))
                
        ), 
        
        tabItem(tabName = "richness",
                
                # richness 
                includeMarkdown("www/richness.md"),
                fluidRow(
                  column(3, uiOutput("yearSelectCombo")),
                  column(3, uiOutput("speciesCombo"))
                ),
                fluidRow(
                  column(6, plotOutput("ggplot2Group3") %>% withSpinner(color = "green")),
                  column(6, plotOutput("ggplot2Group4") %>% withSpinner(color = "green"))
                )
                
        ),
        
        tabItem(tabName = "releases", includeMarkdown("www/releases.md"))
        
      )
      
    ) # end dashboardBody
    
  )# end dashboardPage
  
))