#Load all the necessary packages and libraries.
library(DESeq2)
library(tidyverse)
library(stringr)
library(SummarizedExperiment)
library(purrrlyr)
library(ggplot2)
library(DEGreport)
library(shiny)
library(DT)
library(shinythemes)
shinyUI(fluidPage(
    #Tema de diseño
    theme = shinytheme("united"),
    #Titulo de la aplicación
    titlePanel("mirGFF3 Reader"),
    #Disposición general del interface de la aplicación.
    sidebarLayout(
        #Disposición del panel lateral.
        sidebarPanel(
            #Menú para cargar el arhivo GFF y CSV.
            fileInput("file1", "Choose CSV File",
                      multiple = FALSE,placeholder = "test/metadata.csv",
                      accept = c("text/csv",
                                 "text/comma-separated-values,text/plain",
                                 ".csv")),
            fileInput("file2", "Choose GFF File",
                      multiple = FALSE,placeholder = "test/mirtop.gff",
                      accept = c(".gff")),
            #Casilla activada por defecto que selecciona la normalización o no de los datos
            checkboxInput(inputId = "normalize",
                          label = strong("Data normalization"),
                          value = TRUE),
            #Botón de acción para lanzar la aplicación.
            actionButton("upload", "Upload Data")
        ),
        #Disposición del panel central con un sistema de pestañas.
        mainPanel(
            tabsetPanel(type = "tabs",
                        #Pestaña con información sobre el SummarizedExperiment
                        tabPanel("Info", verbatimTextOutput("contenido")),
                        #Pestaña con la PCA sin colorear
                        tabPanel("PCA", plotOutput("pca"),
                                 #Menu desplegable con las filas de metadata
                                 selectInput("datadrop","metadata", choices = ""),
                                 #Botón de acción para 
                                 actionButton("upload2", "Render plot")
                        ),
                        #Pestaña con una tabla de los datos y la capacidad de seleccionar la fila y obtener los datos de cada isomero.
                        tabPanel("Selector", DT::dataTableOutput("tabla4"),
                                 verbatimTextOutput("info"),plotOutput("graph")),
                        
                        tabPanel("Differential Expression",
                                 hr(),
                                 textInput("formula1", "Enter design", placeholder = "~group"),
                                 p("Type the formula for your differential expression analysis"),
                                 actionButton("upload3", "Do Differential Expression"),
                                 hr(),
                                 selectInput("datadrop2", "Metadata", choices = ""),
                                 actionButton("upload4", "Show results"),
                                 DT::dataTableOutput("tabla5"),
                                 plotOutput("graph2")
                        )
                                
                        
            )
        )
    )
))
