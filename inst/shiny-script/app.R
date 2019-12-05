#' Builds the shiny webapp for DrawAlignR
#'
#' @return None. Calling this script creates the Shiny webapp
#'
#' @import shiny
#' @import ggplot2
#' @import data.table
#' @import plyr
#' @import tibble
#' @import ggplot2
#' @import gridExtra
#' @import ggrepel
#' @import signal
#' @import tidyverse
#' @import crayon
#' @import pbmcapply
#' @import plotly
#' @import mstools
#' @import zoo
#' @import dbplyr
#' @import tidyr
#' @import mzR
#' @import Rcpp
#' @import DIAlignR
#' @import LatticeExtra
#'
#' @source "GetXIC.r"
#' @source "getPepLibData.R"
#' @source "getChromatogramDatapoints.R"
#' @source "plot_aligned.R"
#' @source "plot_chrom_reference.R"

library(plotly)
library(shiny)
library(ggplot2)
library(data.table)
library(plyr)
library(tibble)
library(ggplot2)
library(gridExtra)
library(ggrepel)
library(signal)
library(tidyverse)
library(crayon)
library(pbmcapply)
library(mstools)
library(zoo)
library(dbplyr)
library(tidyr)
library(mzR)
library(Rcpp)
library(DIAlignR)
library(latticeExtra)

#Setting max file size
options(shiny.maxRequestSize=10000*1024^2)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mzR") # Requires netcdf


ui <- fluidPage(

  titlePanel("DrawAlignR"),

  sidebarLayout(

    sidebarPanel(

      #Select 1 or set of mzML files
      fileInput(inputId = "ChromatogramFile", "Choose a Chromatogram File",
                  multiple = TRUE,
                  accept = c(".mzML", ".sqMass")),

      #Select .pqp library file
      fileInput(inputId = "LibraryFile", "Choose a Library File",
                multiple = FALSE,
                accept = c(".PQP")),

      #Path to directory where mzml folder and osw folder are located. By default is set to the working directory.
      textInput(inputId = "WorkingDirectory", "Set Working Directory (Location of mzML and osw folders)", 
                value = paste((gsub('............$', '', getwd())), 'extdata', sep = '')),

      #Full peptide name including modifications
      textInput(inputId = "Mod", "Peptide Name", value = "ANS(UniMod:21)SPTTNIDHLK(UniMod:259)"),

      #Charge of desired peptide (Specific charge must be in data set)
      numericInput(inputId = "Charge", "Peptide Charge", value = 2, min = 1, step = 1),

      #Number of plots to display
      sliderInput("n", "Number of Plots", value=1, min=1, max=10),

      #Off by default. Enabled if DIAlignR should be run and aligned chromatograms should be plotted.
      checkboxInput(inputId = "Align", "Plot Aligned", value = FALSE, width = NULL),

      #Name of the reference run if performing multiple pairwise alignments. Not required.
      textInput(inputId = "Reference", "Select Reference Run for Alignment", value = "chludwig_K150309_013_SW_0"),
    ),

    mainPanel(
      uiOutput("plots")
    )
  )
)

server <- function(input, output) {

  #Generate set of variable plots

  output$plots <- renderUI({
    plot_output_list <- lapply(1:input$n, function(i) {
      plotname <- paste("plot", i, sep="")
      plotlyOutput(plotname)
    })
    do.call(tagList, plot_output_list)
  })

  #Generate all plots. Max plots set to 10

  for (i in 1:10) {
    local({
      my_i <- i
      plotname <- paste("plot", my_i, sep="")

      output[[plotname]] <- renderPlotly({


        #If alignment is disabled, generate standard chromatogram plot.
        
        if (is.null(input$ChromatogramFile)){
          stop()
        }
        else if (is.null(input$LibraryFile)){
          stop()
        }
        else if (is.null(input$Mod)){
          stop()
        }

        if (!(input$Align)){
          chrom_input <- input$ChromatogramFile[[my_i, 'datapath']]
          lib_input <- input$LibraryFile
          peptide <- input$Mod
          lib <- getPepLibData_(lib_input$datapath, peptide_id = '')
          g.out <- getXIC(graphic_obj = ggplot(), chromatogram_file = chrom_input,
                          df_lib = lib, mod = peptide, Isoform_Target_Charge = input$Charge)
          plotly::ggplotly(g.out$graphic_obj, dynamicTicks = TRUE)
        }

        #If alignment is enabled.

        else {

          #Ensuring at least two runs selected, not conducting alignment against same run

          if (!(input$Reference == gsub('...........$', '', input$ChromatogramFile[[my_i, 'name']]))){
            dataPath <- input$WorkingDirectory
            analytes <- paste(input$Mod, "_", toString(input$Charge), sep="")
            runs <- c(input$Reference, gsub('...........$', '', input$ChromatogramFile[my_i, 'name']))
            AlignObjOutput <- getAlignObjs(analytes, runs, dataPath)
            k <- plotAlignedAnalytes(AlignObjOutput, DrawAlignR = TRUE)
            plotly::ggplotly(k[["pBR"]], dynamicTicks = TRUE)
          }

          #If all possible pairwise alignments conducted, plotting the reference run

          else {
            chrom_input <- input$ChromatogramFile[[my_i, 'datapath']]
            lib_input <- input$LibraryFile
            peptide <- input$Mod
            lib <- getPepLibData_(lib_input$datapath, peptide_id = '')
            g.out <- getXIC(graphic_obj = ggplot(), chromatogram_file = chrom_input,
                            df_lib = lib, mod = peptide, Isoform_Target_Charge = input$Charge)
            plotly::ggplotly(g.out$graphic_obj, dynamicTicks = TRUE)
            }
          }
      })
    })
  }
}


shinyApp(ui = ui, server = server)

