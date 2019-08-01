#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

# UI Tabs:
# * Convert
# * Search
# * Process
# * Visualize

library(tidyverse)
library(shiny)
library(shinydashboard)
library(DT)

source('TPP_functions.R')

ui <- dashboardPage(
  dashboardHeader(title="ShinyTPP"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName="home", icon=icon("home", lib="font-awesome")),
      menuItem("Convert", tabName="convert", icon=icon("repeat", lib="glyphicon")),
      menuItem("Search", tabName="search", icon=icon('search', lib="font-awesome")),
      menuItem("Process", tabName="process", icon=icon("filter", lib="font-awesome")),
      menuItem("Visualize & Export", tabName="visualize", icon=icon("sort-by-attributes-alt", lib="glyphicon"))
    )
  ),
  dashboardBody(
    tabItems(
      # Home
      tabItem("home",
              # Tab content goes here
              h2("TPP interface for DDA processing"),
              p("Follow the steps in the sidebar in order to process DDA data with the Trans-proteomic pipeline"),
              p("1. Convert your mass spec raw data into mzML using msconvert"),
              p("2. Search the data using the Comet search engine"),
              p("3. Process the search results to infer proteins and perform label-free quantification"),
              p("4. Visualize the results interactively and export reports"),
              h4("For help, email winget.jm@pg.com")
            ),
      
      # Convert
      tabItem(tabName="convert",
              h2("Convert raw files to mzML"),
              fluidPage(
                box(
              # Tab content goes here
                actionButton("raw.file.select",'Select raw files and convert'),
                dataTableOutput("raw.selected")
                )
              )
            ),
      
      # Search
      tabItem(tabName="search",
              h2("Search data to obtain peptide-spectrum matches"),
              fluidRow(
                # Tab content goes here
                # To generate params we need project dir, db, tol, termini, missed.cleavages,
                # resolution settings, and decoy tag
                box(
                textInput('projectdir', 'Project Directory', 'E:/Project_directory'),
                textInput('database', 'Database', 'E:/dbase/nextprot_all_DECOY.fasta'),
                textInput('decoytag', 'Decoy tag', 'DECOY_'),
                numericInput('tol', 'Precursor tolerance (ppm)', 20),
                numericInput('termini', 'Number of tryptic termini', 1, 0, 2, 1),
                numericInput('missedcleavages', 'Maximum missed cleavages', 2, 0),
                radioButtons('ms2res', 'MS2 resolution', c("High" = "high", "Low" = "low")),
                actionButton('genparams', 'Generate Comet parameters'),
                textOutput('paramsgenerated'),
                title="Generate search parameters"),
                box(
                  actionButton("mzML.file.select", "Select mzML files and search"),
                  dataTableOutput("mzML.selected"),
                  title="Search mzML files"
                )
              )
            ),
    # Process
    tabItem(tabName="process",
            h2("Process search results with the TPP for protein inference and quantification"),
            fluidRow(
              box(
              # Tab content goes here
              textInput('database', 'Database', 'E:/dbase/nextprot_all_DECOY.fasta'),
              checkboxInput('combine', 'Merge into single protXML', value=TRUE),
              textInput('outname', 'Output name if combining', 'combined'),
              actionButton("pepxml.file.select", "Select pepXML files and process")
              )
            )
          ),

  # Visualize
  tabItem(tabName="visualize",
          h2("Examine processed results"),
          fluidRow(
            # Tab content goes here
            actionButton("protxml.file.select", "Select protXML file(s)"),
            dataTableOutput('protxml.extract')
          )
        )
    ) # close tabItems
  ) # close dashboardBody
) # Close dashboardPage

server <- function(input, output, session) {
  # Raw conversion
  observeEvent(input$raw.file.select, {
    raw.files <- choose.files()
    output$raw.selected <- renderDataTable(tibble("Selected files" = raw.files))
    convertRaw(raw.files)
  }, ignoreInit=T, once=T)
  
  # Parameter generation
  observeEvent(input$genparams, {
    if (input$ms2res == 'high'){
      fbt <- 0.02
      fbo <- 0
    } else {
      fbt <- 1.0005
      fbo <- 0.4
    }
    genCometParams(input$projectdir, input$database, input$tol, input$termini, input$missedcleavages, fbt, fbo, input$decoytag)
    output$paramsgenerated <- renderText('Parameters written to disk')
  })
  
  # Search
  observeEvent(input$mzML.file.select, {
    mzml.files <- choose.files()
    output$mzML.selected <- renderDataTable(tibble("Searched files" = mzml.files))
    project.dir = paste0(dirname(mzml.files[1]),'/')
    run.Comet(project.dir, mzml.files)
  }, ignoreInit=T, once=T)
  
  # Process
  observeEvent(input$pepxml.file.select, {
    pepxml.files <- choose.files()
    runTPP(pepxml.files, input$database, input$combine, input$outname)
  }, ignoreInit=T, once=T)
  
  # Visualize
  observeEvent(input$protxml.file.select, {
    protxml.files <- choose.files()
    result <- extract.results(protxml.files)
    output$protxml.extract <- renderDataTable(result, escape=FALSE)
  }, ignoreInit=T, once=T)
}

# Run the application 
shinyApp(ui = ui, server = server)

