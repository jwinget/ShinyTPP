# ShinyTPP

A Shiny app for running the Trans Proteomic Pipeline to search DDA proteomic data

## How to use

1. Clone the repo
2. Open `app.R` in Rstudio
3. Install libraries: tidyverse, shiny, shinydashboard, DT, xml2
4. Click "run App" in the upper-right of app.R window. This should open a web browser
5. Use the tabs on the left for the various TPP functions
  - While they are running, commands will generate output in the Rstudio console so you know they are working
6. Any questions: contact winget.jm@pg.com

## Features to implement

* Multiple search engines. Currently only uses Comet due to an output writing bug in MS-GF+
* "Re-base" two-stage searches for sensitivity with massive databases (e.g. microbiome)
* Additional data visualizations such as descriptive overview graphics and more granular output of processed results
