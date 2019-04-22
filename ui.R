#Paquetes y librerias
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(c("SummarizedExperiment","DESeq2","DEGreport"), version = "3.8")
#library(DESeq2)
#library(tidyverse)
#library(stringr)
#library(SummarizedExperiment)
#library(purrrlyr)
#library(ggplot2)
#library(DEGreport)
library(shiny)

#Iniciamos el Interfaz de usuario para la aplicación shiny.
shinyUI(fluidPage(
  #Título de la aplicación.
  titlePanel("GFF Reader"),
  #Distribución general de la interfaz.
  sidebarLayout(
    #Distribución del panel lateral.
    sidebarPanel(
      #Input de archivos GFF y CSV.
      fileInput("file1", "Choose CSV File"),
      fileInput("file2", "Choose GFF File"),
      #Linea de separación.
      hr(),
      #Selección de columna del archivo de metadatos.
      selectInput("columna", "Columna CSV:",   choices=colnames),
      #Texto explicativo.
      helpText("Choose the CSV column to coloring."),
      #Separación.
      tags$hr(),
      #Cajetilla que registra si se desean normalizar los datos. Desactivada por defecto.
      checkboxInput(inputId = "noise",label = strong("Correcting for noise"),value = FALSE)
    ),
    #Distribución del panel central.
    mainPanel(
      #Diferentes pestañas que muestran las variables "table" y "pca"
      tabsetPanel(type = "tabs",
                  tabPanel("Info", textOutput("info")),
                  tabPanel("Table", plotOutput("table")),
                  tabPanel("PCA", plotOutput("pca"))
      )
    )
  )
))
