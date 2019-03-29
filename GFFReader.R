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
  coldata3<-unlist(str_split(coldata2, ","))
  return(coldata3)
}

##Extracción de COLDATA (opcional) de metadata.csv
# toda la table se necesita para saber que es cada
# muestra.
# La primera columna se toma com0 row.names
# y debe coincidir con `coldata_extract`.
coldata_extract_csv <- function(csv){
  #Cargamos el archivo de metadatos
  md_raw <- read.csv(csv, row.names = 1)
  md_raw
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
  gff_table<-datos1 %>% 
      mutate(uid = str_extract(attribute,"iso-[0-9]+-\\w+;")) %>% # get the UID to have a unique self-identifier
      separate_rows(attribute, sep=";") %>%  # separate by row each element in attribute
      mutate(attribute = trimws(attribute)) %>% # remove leading/tailing spaces
      separate(attribute, sep = " ", into = c("att", "value")) %>% # separate the values into two columns (UID | iso-22-LVMJ3KW9)
      spread(att, value) %>% # move name of the attributes to be columns
      select(-uid) # remove temporal ID
  
  # expresion2<-datos1 %>% separate(attribute, c("Read","UID","Name","Parent","Variant","Cigar","Expression","Filter","Hits"), sep = ";")
  return (gff_table)
}

##Extracción de COUNTS
counts_extract<-function(gff, colnames){
  #Aplicamos la función "atributes_extract" para dividir la variable "Attributes"
  datos1<-atributes_extract(gff)
  #Seleccionamos las variables UID y Expression
  as_tibble(uid_exp <-datos1 %>% 
              select(c(UID, Expression))) %>% 
      separate(Expression, sep=",", into = colnames, convert = TRUE) %>% 
      distinct()  %>%  # no duplicados
      as.data.frame() %>% # because tibble doesn't have row.names
      column_to_rownames("UID") %>%
      as.matrix()
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

gff = "test/mirtop.gff"
metadata = "test/metadata.csv"
colnames<- coldata_extract(gff)
attributes<-atributes_extract(gff)
counts<-counts_extract(gff, colnames)
metadata<-coldata_extract_csv(metadata)
#Creamos la clase SummarizedExperiment
#??
sum_Experiemnt<-SummarizedExperiment(assays = list(raw = counts),
                                     colData = metadata,
                                     rowData = attributes)


# Make function to normalize with DESeq2
# see dds <- DESeq2::DESeqDataSetFromMatrix(counts, metadata, design=~1)
# see vst <- DESeq2::varianceStabilizingTransformation(dds)
# add new normalized count data to new slot in the object 
## sum_Experiemnt@assays[["vst"]] <- assay(vst)

# Make function to PCA
# see DEGreport::degPCA()