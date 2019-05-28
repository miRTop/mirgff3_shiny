
library(shinythemes)
library(DESeq2)
library(tidyverse)
library(stringr)
library(SummarizedExperiment)
library(purrrlyr)
library(ggplot2)
library(DEGreport)
library(shiny)
library(DT)
shinyServer(function(input, output, session) {
    #Función reactiva que procesa el archivo GFF y CSV en un objeto de tipo Summarized Experiment
    options(shiny.maxRequestSize=200*1024^2)
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
                mutate(attribute = gsub("=", " ", trimws(attribute))) %>% # remove leading/tailing spaces
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
        updateSelectInput(session, "datadrop", choices = colnames(metadata))
        keep <- rowSums(counts>0) > (ncol(counts) * 0.2)
        attributes <- attributes[keep,]
        counts <- counts[keep,]
        se1<-SummarizedExperiment(assays = list(raw = counts[,rownames(metadata)]),colData = metadata, rowData = attributes)
        dds<-DESeqDataSetFromMatrix(counts[,rownames(metadata)], metadata,design = ~1)
        vst <- varianceStabilizingTransformation(dds)
        if(input$normalize) {se<-vst} else{se<-se1}
        se
        
    })
    #Condicionamos la salida a que se haya pulsado el boton "upload"
    observeEvent(input$upload2, {
        se<-dataInput()
        output$pca<- renderPlot({
            degPCA(assays(se)[[1]], metadata = colData(se), condition = input$datadrop, data = FALSE)
        })
    })
    
    observeEvent(input$upload, {
        #Definimos la variable "se" (summarizedExperiment) en esta sección para que puedan utilizarla el resto de funciones.
        se<-dataInput()
        #Mostramos el contenido del objeto "contenido" que es un objeto SummarizedExperiment.
        output$contenido<- renderPrint({
           se
        })
        #Mostramos el objeto "pca" que es una PCA del objeto SummarizedExperiment en blanco y negro.
        output$pca<- renderPlot({
            #Blanco y negro
            degPCA(assays(dataInput())[[1]], metadata = colData(se), data = FALSE)
        })
        #Creamos una variable de salida que es una tabla reactiva en la cual podemos seleccionar las filas con la información de rowData y crear gráficos a partir de las lienas seleccionadas.
        #Inicialmente convertimos rowData en un dataFrame para poder hacerlo una tabla ineractiva
        rowdataDF<-as.data.frame(rowData(se))
        #Creamos el output que es una tabla a partir de rowData con lineas seleccionables.
        output$tabla4<-DT::renderDataTable(rowdataDF, server = FALSE )
        #Creamos un output que son gráficos de los isomeros seleccionados.
        output$graph <- renderPlot({
            #Creamos la variable que almacenará las filas seleccionadas.
            filas5 <-input$tabla4_rows_selected
            #Creamos metadata a partir de la variable se
            #Desarrollamos los gráficos de las filas seleccionadas.
            if (!is.null(filas5)){
                degPlot(se, genes = rownames(se)[filas5], slot = 1,
                        xs = input$datadrop, log2 = FALSE)
                
            }
        })
    })
})
