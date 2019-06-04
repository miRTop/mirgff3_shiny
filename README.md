# mirgff3_shiny

a Rshiny tool to visualize mirgff3 files

#Requierements

1) To run the app we need to have R already installed in our computer.


2) Install the package BiocManager that includes:  SummarizedExperiment, DESeq2 and DEGreport

You can type the next code previously to run the app.

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment","DESeq2","DEGreport"))#, version = "3.8")

3) Install also the following packages: tidyverse, shinythemes, stringr, purrrlyr, ggplot2, shiny and DT.

Ej. install.packages(“tidyverse”)

Note* The package "shinythemes" is mandatory for the correct functioning of the program. 
Make sure all the packages are already installed and working in your computer.


#Tests Data


You can find the test data inside mirgff_shiny/test/ 

It consist in a small GFF file (mirtop.gff) and two csv files to choose (metadata.csv, metadata2.csv) to test the app.

#Final Data

In case you want to try the app with a heavier data set you can download the compressed file from:

https://www.dropbox.com/s/addg7wldhl1k46s/geu-mirtop.gff.gz?dl=0

https://www.dropbox.com/s/9xx5d9u1veq6oa3/geu-metadata.csv?dl=0

