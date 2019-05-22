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
shinyUI(fluidPage(
  #Añadimos el tema de la aplicación
  theme = shinytheme("simplex"),
  #Titulo de la aplicación
  titlePanel("mirGFF3 Reader"),
  #Disposición general del interface de la aplicación.
  sidebarLayout(
    #Panel lateral.
    sidebarPanel(
      #Menú para cargar el archivo mirGFF3 y CSV.
      fileInput("file1", "Choose CSV File",
                multiple = FALSE,placeholder = "test/metadata.csv",
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      fileInput("file2", "Choose GFF File",
                multiple = FALSE,placeholder = "test/mirtop.gff",
                accept = c(".gff")),
      
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
                             #Botón de acción para confirmar la elección de la columna de metadata.
                             actionButton("upload2", "Render Plot")
                             ),
                    #Pestaña con una tabla de los datos y la capacidad de seleccionar la fila y obtener los datos de cada isomero.
                    tabPanel("Selector", DT::dataTableOutput("tabla4"),
                             verbatimTextOutput("info"),plotOutput("graph")),
                    #Pestaña con información sobre expresión diferencial
                    tabPanel("Differential Expression", 
                             #Casilla activada por defecto que selecciona la normalización o no de los datos
                             checkboxInput(inputId = "normalize",
                                           label = strong("Data normalization"),
                                           value = TRUE),
                             #Encabezado de DESeqDataSet normalizado (varianceStabilizingTransformation) o no.
                             verbatimTextOutput("expresion"), 
                             hr(),
                             #Botón que permite mostrar los valores de dispersion.
                             actionButton("upload3", "Dispersion"),
                             verbatimTextOutput("expresion2"),
                             hr(),
                             #Boton que permite mostrar el DEGplot de los 12 primeros isomiR.
                             actionButton("upload4", "Grafico"),
                             plotOutput("expresion_plot"),
                             #Boton que permite mostrar el MA-Plot.
                             actionButton("upload5", "Grafico2"),
                             plotOutput("expresion_plot2")
                             )

        )
    )
  )
))
