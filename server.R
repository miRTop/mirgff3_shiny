
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
shinyServer(function(input, output) {
  #Función reactiva que procesa el archivo GFF y CSV en un objeto de tipo Summarized Experiment
  dataInput <-reactive({
    inFile1 <-input$file1
    inFile2 <-input$file2
    lecturaGFF <- function(gff){
      resultado<- as_tibble(read_tsv(gff, col_names = c("seqname","source","feature","start","end","score","strand","frame","attribute"), comment = "#"))
      return (resultado)
    } 
    coldata_extract <- function(gff_coldata) {
      primeras_lineas = readLines(gff_coldata, n = 3)
      coldata <- grep("COLDATA",primeras_lineas, value = TRUE)
      coldata2<-coldata %>% str_replace("## COLDATA: ","")
      coldata3<-unlist(str_split(coldata2, ","))
      return(coldata3)
    }
    coldata_extract_csv <- function(csv){
      md_raw <- read.csv(csv, row.names = 1)
      return(md_raw)
    }
    atributes_extract<-function(gff){
      datos1<-lecturaGFF(gff)
      as_tibble(expresion <-datos1 %>% 
                  select(attribute))
      gff_table<-datos1 %>% 
        mutate(uid = str_extract(attribute,"iso-[0-9]+-\\w+;")) %>% # get the UID to have a unique self-identifier
        separate_rows(attribute, sep=";") %>%  # separate by row each element in attribute
        mutate(attribute = trimws(attribute)) %>% # remove leading/tailing spaces
        separate(attribute, sep = " ", into = c("att", "value")) %>% # separate the values into two columns (UID | iso-22-LVMJ3KW9)
        spread(att, value) %>% # move name of the attributes to be columns
        select(-uid) # remove temporal ID
      return (gff_table)
    }
    counts_extract<-function(gff, colnames){
      datos1<-atributes_extract(gff)
      as_tibble(uid_exp <-datos1 %>% 
                  select(c(UID,Expression))) %>% 
        select(c(UID, Expression)) %>% 
        separate(Expression, sep=",", into = colnames, convert = TRUE) %>% 
        distinct()  %>%  # no duplicados
        as.data.frame() %>% # because tibble doesn't have row.names
        column_to_rownames("UID") %>%
        as.matrix()
    }
    colnames<- coldata_extract(inFile2$datapath)
    attributes<-atributes_extract(inFile2$datapath)
    counts<-counts_extract(inFile2$datapath, colnames)
    metadata<-coldata_extract_csv(inFile1$datapath)
    se<-SummarizedExperiment(assays = list(raw = counts), colData = metadata, rowData = attributes)
    se
  })
  #Condicionamos la salida a que se haya pulsado el boton "upload"
  observeEvent(input$upload, {
    #Mostramos el contenido del objeto "contenido" que es un objeto SummarizedExperiment.
    output$contenido<- renderPrint({
      dataInput()
    })
    #Mostramos el objeto "pca" que es una PCA del objeto SummarizedExperiment en blanco y negro.
    output$pca<- renderPlot({
      degPCA(assays(se)[[1]], metadata = colData(se), data = FALSE)
    })
  })
})
