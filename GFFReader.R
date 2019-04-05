#Paquetes y librerias
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment","DESeq2"), version = "3.8")
library(DESeq2)
library(tidyverse)
library(stringr)
library(SummarizedExperiment)
library(purrrlyr)


###EXTRACCION DE DATOS

####EXTRACCION DE DATOS

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
##Extracción de la columna ATRIBUTOS/ROWDATA
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
  
  # expresion2<-datos1 %>% separate(attribute, c("Read","UID","Name","Parent","Variant","Cigar","Expression","Filter","Hits"), sep = ";")
  return (gff_table)
}

##Extracción de COUNTS
  counts_extract<-function(gff, colnames){
    #Aplicamos la función "atributes_extract" para dividir la variable "Attributes"
    datos1<-atributes_extract(gff)
    #Seleccionamos las variables UID y Expression
    as_tibble(uid_exp <-datos1 %>% 
                select(c(UID,Expression)))
    #Transformamos el dibble en matriz
    uid_exp2<-as.matrix(uid_exp)
    return(uid_exp)
    select(c(UID, Expression)) %>% 
  separate(Expression, sep=",", into = colnames, convert = TRUE) %>% 
  distinct()  %>%  # no duplicados
  as.data.frame() %>% # because tibble doesn't have row.names
  column_to_rownames("UID") %>%
  as.matrix()
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


##PRUEBA
#Cargamos los archivos necesarios
# coldataGFF = ERR187525-mirbase-ready ERR187844-mirbase-ready ERR187497-mirbase-ready ERR187773-mirbase-ready
coldataGFF<-coldata_extract("test/mirtop.gff")
# coldataGFF2 = ERR187525-mirbase-ready ERR187844-mirbase-ready ERR187497-mirbase-ready ERR187773-mirbase-ready
coldataGFF2<-coldata_extract_csv("test/metadata.csv")
#rowdata1 = seqname | source | feature | start | end | score | strand | frame | ATTRIBUTE
#                                                 Read | UID | Name | Parent | Variant | Cigar | Expression | Filter | Hits                        
rowdata1<-atributes_extract("test/mirtop.gff")
#counts1 = UID - Expression
counts1<-counts_extract("test/mirtop.gff",coldataGFF)

#Aplicación de SummarizedExperiment
sum_Experiment<-SummarizedExperiment(assays = list(raw = counts1),colData = coldataGFF2,rowData = rowdata1)

#Error in all_dims[, 1L] : número incorreto de dimensiones

dds<- DESeqDataSetFromMatrix(countData = counts1, colData = coldataGFF2, design = ~1)
#He intentado aplicar esta función pero no acabo de ver como hacerlo. Es necesatrio que el numero de columnas de count1 (UID, Expression) sean igual al numero
#de lineas de coldataGFF2 (ERR187525-mirbase-ready, ERR187844-mirbase-ready, ERR187497-mirbase-ready, ERR187773-mirbase-ready).
#No se si te refieres al numero en Expression con respecto a los valores de coldataGFF2
#Por si es eso:
#Eliminamos el valor UID
counts2<-counts1[,-1]
#Lo convertimos en un data frame para trabajar con el.
counts3 <- as.data.frame(counts2)
#Dividimos los valores de la columna de Expression y creamos una columna por cada gupo de valores
valor <- strsplit(as.character(counts3$Expression), ",")
counts3$col1 <- sapply(valor, `[`, 1)
counts3$col2 <- sapply(valor, `[`, 2)
counts3$col3 <- sapply(valor, `[`, 3)
counts3$col4 <- sapply(valor, `[`, 4)
#Eliminamos la columna de expression
counts4<-counts3[,-1]
#Cambiamos los nombres de las columnas por los valores de metadata
colnames(counts4)<- coldataGFF
#Cambiamos el tipo de dato de character a numerico
counts4[] <- lapply(counts4, function(x) as.numeric(as.character(x)))
#Aplicamos el DESeqDataSetFromMatrix
dds<- DESeqDataSetFromMatrix(countData = counts4, colData = coldataGFF2, design = ~1)
vst <- DESeq2::varianceStabilizingTransformation(dds)
