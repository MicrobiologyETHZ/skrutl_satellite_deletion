library(DESeq2)
library(tximport)
library(edgeR)
library(tidyverse)
source("code/deseq_fx.R")


sampleData_file <- "data/deletion/22-02-23-deletion-sample-data.csv"
conditions <- "Genotype"
#tx2gene_file <- "data/salmon/03_22/22-08-22-Drosophila_melanogaster.BDGP6.32.104.tx2gene.csv"
outDir <- "data/deletion/results"
prefix <- "22-02-23-deletion-featCnts-"
lfct <- 0.5
alpha <- 0.01
annotations <- "data/22-08-22-org.Dm.eg.db-annotation.csv"
count_file <- "data/deletion/BDGP6_deletion.merged.featureCounts.csv"

deseq_output <- run_deseq_on_featCnts(sampleData_file, count_file, conditions,
                                      outDir, prefix)

dds <- deseq_output$dds
sd <- deseq_output$sampleData


pattern1 <- '^M|Z'
pattern2 <- 'yw'
get_all_results(sd, pattern1, pattern2, dds, outDir, prefix, annotations)



pattern1 <- 'Zhr_M'
pattern2 <- '^M'

get_all_results(sd, pattern1, pattern2, dds, outDir, prefix, annotations)


pattern1 <- '^M|Zhr_'
pattern2 <- 'Zhr'

get_all_results(sd, pattern1, pattern2, dds, outDir, prefix, annotations)

pattern1 <- 'Zhr_M41A10'
pattern2 <- 'Zhr_NipD'

get_all_results(sd, pattern1, pattern2, dds, outDir, prefix, annotations)

