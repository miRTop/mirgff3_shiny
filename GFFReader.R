#Paquetes y librerias
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment","DESeq2"), version = "3.8")
library(DESeq2)
library(tidyverse)
library(stringr)
library(SummarizedExperiment)
library(purrrlyr)

####EXTRACCION DE DATOS DESDE FICHEROS GFF Y CSV

##Función general de lectura de GFF
#seqname | source | feature | start | end | score | strand | frame | attribute
lecturaGFF <- function(gff){
  resultado<- as_tibble(read_tsv(gff, col_names = c("seqname","source","feature","start","end","score","strand","frame","attribute"), comment = "#"))
  return (resultado)
} 
##Extracción de ## COLDATA
#ERR187525-mirbase-ready ERR187844-mirbase-ready ERR187497-mirbase-ready ERR187773-mirbase-ready
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
coldata_extract_csv <- function(csv){
  #Cargamos el archivo de metadatos
  md_raw <- read.csv(csv, row.names = 1)
  return(md_raw)
}

###PROCESAMIENTO DE DATOS A PARTIR DE LO EXTRAIDO DE LOS FICHEROS

##Extracción de la columna TTRIBUTE y creación de ROWDATA
atributes_extract<-function(gff){
  #Leemos el archivo GFF
  datos1<-lecturaGFF(gff)
  #Seleccionamos la columna "attribute" y almaceamos sus valores en datos1
  as_tibble(expresion <-datos1 %>% 
              select(attribute))
  #Modificamos los atributos que constituyen datos1
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
    #Aplicamos la función "atributes_extract" para dividir la variable "Attributes"
    datos1<-atributes_extract(gff)
    #Seleccionamos las variables UID y Expression
    as_tibble(uid_exp <-datos1 %>% 
                select(c(UID,Expression))) %>% 
    #Transformamos el dibble en matriz
    # uid_exp2<-as.matrix(uid_exp)
    # return(uid_exp)
        select(c(UID, Expression)) %>% 
        separate(Expression, sep=",", into = colnames, convert = TRUE) %>% 
        distinct()  %>%  # no duplicados
        as.data.frame() %>% # because tibble doesn't have row.names
        column_to_rownames("UID") %>%
        as.matrix()
  }

###ASIGNACIÓN Y CREACION DE VARIABLES
gff = "test/mirtop.gff"
metadata = "test/metadata.csv"
colnames<- coldata_extract(gff)
attributes<-atributes_extract(gff)
counts<-counts_extract(gff, colnames)
metadata<-coldata_extract_csv(metadata)


#Aplicación de SummarizedExperiment
se<-SummarizedExperiment(assays = list(raw = counts),
                         colData = metadata,
                         rowData = attributes)

#Alamcenamos ls valore siniciales y finales del análisis de expresión diferencial.
dds<- DESeqDataSetFromMatrix(counts, metadata, design = ~1)

#Calculamos los valores de dispersión.
dds <- estimateSizeFactors(dds)

#Normalización de los parámetros.
vst <- varianceStabilizingTransformation(dds)

#Creamos el assay "vst" dentro del objeto SummarizedExperiment.
assays(se)[["vst"]] = assay(vst)

#Dibujamos el plot de los 12 primeros isomeros de datos crudos.
degPlot(dds,genes = rownames(dds)[1:12], xs = "group",log2 = FALSE)

#Dibujamos el plot de los 12 primeros isomeros de datos normalizados.
degPlot(vst,genes = rownames(vst)[1:12], xs = "group",log2 = FALSE)

#Realizamos la PCA coloreando las observacionen según al grupo que pertenecen.
degPCA(counts, metadata = metadata, condition = "group", data = FALSE)


# Make function to PCA
# http://lpantano.github.io/DEGreport/reference/degPCA.html
# see DEGreport::degPCA(se) #PCA

# http://lpantano.github.io/DEGreport/reference/degPlot.html
## xs = column of colData(se)
## slot = "vst"
## log2 = FALSE
# see DEGreport::degPlot(se) # plot of specific isomirs
