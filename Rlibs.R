args = commandArgs(trailingOnly=TRUE)
libloc <-  args[1]
reposloc<-"https://cloud.r-project.org/" #consider args[2]
dir.create(libloc)

.libPaths(libloc)
if (!require("devtools")) install.packages("devtools", lib = libloc, repos = reposloc)
if (!require("remotes")) install.packages("remotes", lib = libloc, repos = reposloc)
if (!require("withr")) install.packages("withr", lib = libloc, repos = reposloc)
if (!require("vctrs")) install.packages("vctrs", lib = libloc, repos = reposloc)
if (!require("patchwork")) with_libpaths(new = libloc, devtools::install_github("thomasp85/patchwork"))
if (!require("dsb")) install.packages("dsb", lib = libloc, repos = reposloc) #to verify if same
if (!require("tidyverse")) install.packages("tidyverse", lib = libloc, repos = reposloc)
if (!require("magrittr")) install.packages("magrittr", lib = libloc, repos = reposloc)
if (!require("data.table")) install.packages("data.table", lib = libloc, repos = reposloc)
if (!require("Matrix")) install.packages("Matrix", lib = libloc, repos = reposloc)
if (!require("ggplot2")) install.packages("ggplot2", lib = libloc, repos = reposloc)
if (!require("readr")) install.packages("readr", lib = libloc, repos = reposloc)
if (!require("Seurat")) install.packages("Seurat", lib = libloc, repos = reposloc)
if (!require("intrinsicDimension")) install.packages("intrinsicDimension", lib = libloc, repos = reposloc)
if (!require("scater")) install.packages("scater", lib = libloc, repos = reposloc)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",lib = libloc, repos = reposloc)
if (!require("DropletUtils")) BiocManager::install("DropletUtils", lib = libloc)
if (!require("scry")) BiocManager::install("scry", lib = libloc) #requires R >=4.0.3
