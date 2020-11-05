#!/usr/bin/env Rscript

library(shiny)

#' Launch EPItOMe application.
#'
#' Launches the Shiny web interface for EPItOMe
#'
#' @import shinydashboard
#' @import DT
#' @import R.utils
#' @importFrom shinyjs enable disable toggle info onclick hide show
#' @importFrom shinyBS addPopover
#' @importFrom shiny runApp
#' @importFrom Biostrings readAAStringSet AAStringSet
#' @importFrom utils data
#' @export
#' 

#install.packages("hex", repos = NULL, type="source")
#launch_EPItOMeSim <- function(){
#	appDir <- system.file("EPItOMe-app", package = "epitomeSim")
#	print(system.file("EPItOMe-app", package = "epitomeSim"))#
#	if (appDir == "") {
#	stop("Unable to launch INfORM Shiny interface, please check installation!.", call. = FALSE)
#	}
#
#	shiny::runApp(appDir, display.mode = "normal")
#}
#launch_EPItOMeSim()
