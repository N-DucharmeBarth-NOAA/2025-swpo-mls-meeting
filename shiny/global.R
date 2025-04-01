# global.R - Load packages and data that are used across UI and server

# Load required packages
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(data.table)
library(markdown)
library(magrittr)
library(ggplot2)
library(viridis)
library(ggthemes)
library(fresh)

# Load data once to improve performance
  summary_dt = fread(file.path("data","summary_dt.csv"))

# show available directories
    all_dirs = list.files("data",recursive = TRUE)
    all_dirs = all_dirs[grep("html_parameters.txt",all_dirs,fixed=TRUE)]
    all_dirs = gsub("html_parameters.txt","",all_dirs,fixed=TRUE)

    model_stem = "data"
