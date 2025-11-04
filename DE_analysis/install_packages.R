# Install BiocManager if needed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c("DESeq2", "tximport", "edgeR"))

# Install tidyverse from CRAN
install.packages("tidyverse")