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

#Iniciamos el archivo "server" de la aplicación shiny.
shinyServer(function(input, output) {
  
##########Código para porcesar archivos GFF##########
  
  #Leer archivo GFF
  lecturaGFF <- function(gff){
    resultado<- as_tibble(read_tsv(gff, col_names = c("seqname","source","feature","start","end","score","strand","frame","attribute"), comment = "#"))
    return (resultado)
  } 
  coldata_extract <- function(gff_coldata) {
    #Leemos las lineas de comentarios.
    primeras_lineas = readLines(gff_coldata, n = 3)
     #Capturamos la linea con el texto "COLDATA"
   coldata <- grep("COLDATA",primeras_lineas, value = TRUE)
    #Eliminamos el texto superfluo
    coldata2<-coldata %>% str_replace("## COLDATA: ","")
    #Separamos los distintos elementos
    coldata3<-unlist(str_split(coldata2, ","))
  return(coldata3)
}
  #Extracción de COLDATA
  coldata_extract_csv <- function(csv){
    md_raw <- read.csv(csv, row.names = 1)
    return(md_raw)
  }
  #Extracción de ROWDATA
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
  ##Extracción de COUNTS
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
  ##########Fin del bloque de código para procesar GFF##########
  
  #Función PROVISIONAL
  #Texto con los datos del archivo procesado GFF y CSV. Sale como "output$info"
  output$info<- renderText({
    #Almacenamos el archivo CSV como inFile1
    inFile1 <-input$file1
    if(is.null(inFile1))
      return(NULL)
    #Alacenamos el archiv GFF como inFile2
    inFile2 <-input$file2
    if(is.null(inFile2))
      return(NULL)
    #Establecemos el valor de: rowdata, coldata, metadata y counts.
    colnames<- coldata_extract(inFile2$datapath)
    attributes<-atributes_extract(inFile2$datapath)
    counts<-counts_extract(inFile2$datapath, colnames)
    metadata<-coldata_extract_csv(inFile1$datapath)
    #Creamos el objeto SumamarizedExperiment
    se<-SummarizedExperiment(assays = list(raw = counts), colData = metadata, rowData = attributes)
    #Mostramos la información en el Assay "raw"
    head(assays(se)[["raw"]])
    
    #En caso de que se marque la casilla de normalización.
    if (input$noise) {
      #Normalizamos los datos.
      dds<- DESeqDataSetFromMatrix(counts, metadata, design = ~1)
      vst<-varianceStabilizingTransformation(dds)
      #Creamos el Assay "vst" dentro del objeto SummarizedExperiment
      assays(se)[["vst"]] = assay(vst)
      #Mostramos la información en el Assay "vst"
      head(assays(se)[["vst"]])
      }
  })
  
  #Creamos la variable de salida table que son plots con la información de GFF y CSV.
  output$table<- renderPlot({
    #Almacenamos el archivo CSV como inFile1
    inFile1 <-input$file1
    if(is.null(inFile1))
      return(NULL)
    #Alacenamos el archiv GFF como inFile2
    inFile2 <-input$file2
    if(is.null(inFile2))
      return(NULL)
    #Establecemos el valor de: rowdata, coldata, metadata y counts.
    colnames<- coldata_extract(inFile2$datapath)
    attributes<-atributes_extract(inFile2$datapath)
    counts<-counts_extract(inFile2$datapath, colnames)
    metadata<-coldata_extract_csv(inFile1$datapath)
    #Registramos la columna de CSV seleccionada como "seleccion"
    seleccion<-input$columna
    #Creamos el objeto SumamarizedExperiment
    se<-SummarizedExperiment(assays = list(raw = counts), colData = metadata, rowData = attributes)
    dds<- DESeqDataSetFromMatrix(counts, metadata, design = ~1)
    dds <- estimateSizeFactors(dds)
    degPlot(dds,genes = rownames(dds)[1:12], xs = "group",log2 = FALSE)
    
    #En caso de que se marque la casilla de normalización.
    if (input$noise) {
      #Normalizamos los datos.
      dds<- DESeqDataSetFromMatrix(counts, metadata, design = ~1)
      vst<-varianceStabilizingTransformation(dds)
      #Creamos el Assay vst dentro del objeto SummarizedExperiment.
      #assays(se)[["vst"]] = assay(vst)
      degPlot(vst,genes = rownames(vst)[1:12], xs = "group",log2 = FALSE)
    }
  })
  output$pca<- renderPlot({
    degPCA(counts, metadata = metadata, condition = "group", data = FALSE)
  })
})
