args = commandArgs(trailingOnly=TRUE)
libloc <-  args[1]

dir.create(libloc)

.libPaths(libloc)
if (!require("devtools")) install.packages("devtools", lib = libloc)
if (!require("withr")) install.packages("withr", lib = libloc)
if (!require("vctrs")) install.packages("vctrs", lib = libloc)
if (!require("patchwork")) with_libpaths(new = libloc, install_github("thomasp85/patchwork"))
if (!require("dsb")) with_libpaths(new = libloc, install_github("MattPM/dsb"))
if (!require("tidyverse")) install.packages("tidyverse", lib = libloc)
if (!require("magrittr")) install.packages("magrittr", lib = libloc)
if (!require("data.table")) install.packages("data.table", lib = libloc)
if (!require("Matrix")) install.packages("Matrix", lib = libloc)
if (!require("ggplot2")) install.packages("ggplot2", lib = libloc)
if (!require("readr")) install.packages("readr", lib = libloc)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager",lib = libloc)
if (!require("DropletUtils")) BiocManager::install("DropletUtils", lib = libloc)

