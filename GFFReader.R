#Paquetes y librerias
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment", version = "3.8")
library(tidyverse)
library(stringr)
library(SummarizedExperiment)

###EXTRACCION DE DATOS

##Función general de lectura de gff en seqname | source | feature | start | end | score | strand | frame | attribute
lecturaGFF <- function(gff){
  resultado<- as_tibble(read_tsv(gff, col_names = c("seqname","source","feature","start","end","score","strand","frame","attribute"), comment = "#"))
  return (resultado)
} 

##Extracción de ## COLDATA
coldata_extract <- function(gff_coldata) {
  #Leemos las lineas de comentarios.
  primeras_lineas = readLines(gff_coldata, n = 3)
  #Capturamos la linea con el texto "COLDATA"
  coldata <- grep("COLDATA",primeras_lineas, value = TRUE)
  #Eliminamos el texto superfluo
  coldata2<-coldata %>% str_replace("## COLDATA: ","")
  #Separamos los distintos elementos
  coldata3<-str_split(coldata2, ",")
  return(coldata3)
}

##Extracción de COLDATA (opcional) de metadata.csv
coldata_extract_csv <- function(csv){
  #Cargamos el archivo de metadatos
  md_raw <- read_csv(csv)
  #Seleccionamos los samplenames
  samplename <- md_raw %>% pull(samplename)
  return(samplename)
}

##Extracción ROWDATA (GFF sin la columna atribute )
rowdata_extract<- function(gff){
  #Empleamos la función de lectura de GFF
  resultado<-lecturaGFF(gff)
  #Tomamos todas las columnas excepto "Attribute"
  as_tibble(resto <-resultado %>% 
              select(-attribute))
  return(resto)
}

##Extracción de la columna ATRIBUTOS
atributes_extract<-function(gff){
  #Leemos el archivo GFF
  datos1<-lecturaGFF(gff)
  #Seleccionamos la columna "attribute" y almaceamos sus valores
  as_tibble(expresion <-datos1 %>% 
              select(attribute))
  #Separamos los distintos valores de la linea de "attribute" en sus diferentes tipos 
  expresion2<-datos1 %>% separate(attribute, c("Read","UID","Name","Parent","Variant","Cigar","Expression","Filter","Hits"))
  return (expresion2)
}

##Extracción de COUNTS
counts_extract<-function(gff){
  #Aplicamos la función "atributes_extract" para dividir la variable "Attributes"
  datos1<-atributes_extract(gff)
  #Seleccionamos las variables UID y Expression
  as_tibble(uid_exp <-datos1 %>% 
              select(c(UID,Expression)))
  #Transformamos el dibble en matriz
  uid_exp2<-as.matrix(uid_exp)
  return(uid_exp)
}

##Función Englobadora. Introducimos un GFF y devuelve un objeto con tres bloques de información: counts, rowdata y coldata

gff_objects<-function(gff){
  #Aplicamos las diversas funciones para obtener los tres objetos
  counts1<-counts_extract(gff)
  rowdata1<-rowdata_extract(gff)
  coldata1<- coldata_extract(gff)
  #Creamos una lista con los resultados
  newList <- list(counts1, rowdata1, coldata1)
  return(newList)
}

#Creamos la clase SummarizedExperiment
#??
sum_Experiemnt1<-makeSummarizedExperimentFromExpressionSet(gff_objects())
