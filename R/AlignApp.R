options(shiny.maxRequestSize=10000*1024^2)

library(shiny)
library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(signal)
library(tidyverse)

ui <- fluidPage(
  titlePanel("DrawAlignR"),
  sidebarLayout(
    sidebarPanel(
      #sliderInput("RT_Range", "RT_Range", 0, 5000, c(25, 40)),
      fileInput(inputId = "ChromatogramFile", "Choose a Chromatogram File",
                  multiple = TRUE,
                  accept = c(".mzML", ".sqMass")),
      fileInput(inputId = "LibraryFile", "Choose a Library File",
                multiple = TRUE,
                accept = c(".PQP")),
      textInput(inputId = "Mod", "Peptide Name", value = "")
    ),
    mainPanel(
      plotOutput("Chromatogram")
    )
  )
)

server <- function(input, output) {
  output$Chromatogram <- renderPlot({
    lib <- getPepLibData_(input$LibraryFile, peptide_id = '')
    g <- ggplot()
    g.out <- getXIC(graphic_obj = g, chromatogram_file = input$ChromatogramFile, df_lib = lib, mod =input$Mod)
  })
}

shinyApp(ui = ui, server = server)
