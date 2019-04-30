#Load all the necessary packages and libraries.
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment","DESeq2","DEGreport"), version = "3.8")
library(DESeq2)
library(tidyverse)
library(stringr)
library(SummarizedExperiment)
library(purrrlyr)
library(ggplot2)
library(DEGreport)
library(shiny)
shinyUI(fluidPage(
  titlePanel("mirGFF3 Reader"),
  sidebarLayout(
    sidebarPanel(
      #Menú para cargar el arhivo GFF y CSV.
      fileInput("file1", "Choose CSV File",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      fileInput("file2", "Choose GFF File",
                multiple = FALSE,
                accept = c(".gff")),
      #Botón de acción para lanzar la aplicación.
      actionButton("upload", "Upload Data")
    ),
    #Panel central con un sistema de pestañas.
    mainPanel(
        tabsetPanel(type = "tabs",
                    #Pestaña con información sobre el SummarizedExperiment
                    tabPanel("Info", verbatimTextOutput("contenido")),
                    #Pestaña con la PCA sin colorear
                    tabPanel("PCA", plotOutput("pca"))
        )
    )
  )
))
