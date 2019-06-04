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
