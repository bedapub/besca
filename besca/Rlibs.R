args = commandArgs(trailingOnly=TRUE)
libloc <-  args[1]

dir.create(libloc)

.libPaths(libloc)
if (!require("devtools")) install.packages("devtools", lib = libloc)
if (!require("withr")) install.packages("withr", lib = libloc)
if (!require("vctrs")) install.packages("vctrs", lib = libloc)
if (!require("patchwork")) install.packages("patchwork", lib = libloc)
if (!require("dsb")) with_libpaths(new = libloc, install_github("MattPM/dsb"))
if (!require("tidyverse")) install.packages("tidyverse", lib = libloc)
if (!require("magrittr")) install.packages("magrittr", lib = libloc)
if (!require("data.table")) install.packages("data.table", lib = libloc)
if (!require("Matrix")) install.packages("Matrix", lib = libloc)
if (!require("DropletUtils")) with_libpaths(new = libloc, install_github("MarioniLab/DropletUtils"))
if (!require("BiocManager")) install.packages("BiocManager", lib = libloc)
if (!require("scater")) BiocManager::install("scater", lib = libloc)

