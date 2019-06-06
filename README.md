# mirgff3_shiny

a Rshiny tool to visualize mirgff3 files

## dependencies

```
pkgs <- c(
    "shinythemes",
    "DESeq2",
    "tidyverse",
    "stringr",
    "SummarizedExperiment",
    "purrrlyr",
    "ggplot2",
    "DEGreport",
    "shiny",
    "DT")
    
install.packages("BiocManager)
BiocManager::install(pkgs_to_install, update=FALSE, ask=FALSE)
```

## Quick start

  * [Download app](https://github.com/miRTop/mirgff3_shiny/archive/master.zip)
  * Extract file
  * Open Rproj file

```
library(shinythemes)
runApp()
```

Use the files in test to use the app with a dummy data.

## Optional data

For a test with heavier data files download the compressed files from:

https://www.dropbox.com/s/addg7wldhl1k46s/geu-mirtop.gff.gz?dl=0

https://www.dropbox.com/s/9xx5d9u1veq6oa3/geu-metadata.csv?dl=0


