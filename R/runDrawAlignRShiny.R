#' Launches the shiny app for DrawAlignR - User facing function.
#'
#'
#' @return None. Runs the Shiny webapp found in inst/available-shiny-apps/DrawAlignR/app.R
#'
#'
#'
#' @export
#' @import shiny


runDrawAlignR <- function() {
  appDir <- system.file("ShinyScript",
                        package = "DrawAlignR")

  shiny::runApp(appDir, display.mode = "normal")
  return()
}


