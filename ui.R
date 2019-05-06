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
  titlePanel("mirGFF3 Reader"),
  sidebarLayout(
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
      #Botón de acción para lanzar la aplicación.
      actionButton("upload", "Upload Data")
    ),
    #Panel central con un sistema de pestañas.
    mainPanel(
        tabsetPanel(type = "tabs",
                    #Pestaña con información sobre el SummarizedExperiment
                    tabPanel("Info", verbatimTextOutput("contenido")),
                    #Pestaña con la PCA sin colorear
                    tabPanel("PCA", plotOutput("pca")),
                    #Pestaña con la tabla de isomeros y un gráfico con los mismos donde podemos seleccionar las filas a destacar.
                    tabPanel("Prueba", fluidRow(column(6, DT::dataTableOutput("tabla1")),
                                               column(6, plotOutput("grafico1"))
                                                )
                    ),
                    #Pestaña con la tabla de ismoeros y un gráfico con los mismos donde se destacan los que aparecen en la tabla actual.
                    tabPanel("Prueba2", fluidRow(column(6, DT::dataTableOutput("tabla2")),
                                                column(6, plotOutput("grafico2"))
                                                )
                    ),
                    #Pestaña con la tabla de isomeros y un gráficos de los mismos donde podemos seleccionar los que aparecen en pantalla por columnas.
                    tabPanel("Prueba3", fluidRow(column(6, DT::dataTableOutput("tabla3")),
                                                 column(6, plotOutput("grafico3"))
                                                )
                    ),
                    #Pestaña con una tabla de los datos y la capacidad de seleccionar la fila y obtener los datos de cada isomero.
                    tabPanel("Selector", DT::dataTableOutput("tabla4"),
                             verbatimTextOutput("info"))

        )
    )
  )
))
