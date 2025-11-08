########################### load general dependencies ########################### 
packagesS<-c("readxl",
             "dplyr",
             "magrittr",
             "tibble",
             "limma",
             "DESeq2",
             "edgeR",
             "ggplot2",
             "forcats",
             "tidyr",
             "RColorBrewer",
             "ggvenn",
             "gplots",
             "grDevices",
             "viridis",
             "ComplexHeatmap",
             "tibble",
             "cluster",
             "stringr",
             "factoextra",
             "circlize",
             "extrafont"
             
             
             
             
)

invisible(suppressWarnings(suppressMessages(package.check <- lapply(
  packagesS,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
))))

font_import()