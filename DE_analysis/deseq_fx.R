#data_dir: salmon file directory
#library(DESeq2)


run_deseq_on_featCnts <- function(sampleData_file, count_file, 
                                  conditions, outDir, prefix, rundeseq=TRUE){
  
  # Validating the data
  
  sampleData <- read.csv(sampleData_file, header=TRUE, row.names=1)
  samples <- rownames(sampleData)
  conditions <- unlist(strsplit(conditions, ','))
  sampleData <- sampleData %>% unite('group', conditions, remove=FALSE)
  
  # Importing the data
  cnts <- read_csv(count_file)
  cnts <- cnts %>% select(-Length) %>% as.data.frame()
  rownames(cnts) <- cnts$Geneid
  cnts = cnts[, samples]
  # Differential analysis
  
  if (rundeseq == TRUE){
    dds <- DESeqDataSetFromMatrix(countData = cnts,
                                  colData = sampleData,
                                  design = ~ group)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    dds <- estimateSizeFactors(dds)
    norm_cnts <- counts(dds, normalized=TRUE)
    dds <-  DESeq(dds)
    vsd <- vst(dds)
    write.csv(assay(vsd), file.path(outDir, paste0(prefix, "featCnts-vsd.csv")))
    
    write.csv(norm_cnts, file.path(outDir, paste0(prefix, "norm_cnts.csv")))
    
    return(list("dds"=dds, "sampleData"=sampleData))
    
  }
  return(list())
}





run_deseq_on_salmon <- function(sampleData_file, data_dir, tx2gene_file,
                      conditions, outDir, prefix, rundeseq=TRUE){
  
  # Validating the data
  
  sampleData <- read.csv(sampleData_file, header=TRUE, row.names=1)
  samples <- rownames(sampleData)
  sample_files <- file.path(data_dir,  paste0(samples, "_quant"), "quant.sf")
  names(sample_files) <- samples
  tx2gene <- read.csv(tx2gene_file)
  conditions <- unlist(strsplit(conditions, ','))
  
  sampleData <- sampleData %>% unite('group', conditions, remove=FALSE)
  
  # Importing the data
  txi <- tximport(sample_files, type = "salmon", tx2gene = tx2gene)
  print(dim(txi$counts))
  write.csv(txi$abundance, file.path(outDir, paste0(prefix, "salmon-gene-tpms.csv")))
  
  # Differential analysis
  
  if (rundeseq == TRUE){
    dds <- DESeqDataSetFromTximport(txi, sampleData, ~ group)
    keep <- rowSums(counts(dds)) >= 10
    dds <- dds[keep,]
    dds <- DESeq(dds)
    vsd <- vst(dds)
    
    write.csv(assay(vsd), file.path(outDir, paste0(prefix, "salmon-vsd.csv")))
    return(list("dds"=dds, "sampleData"=sampleData))
    
  }
  return(list())
}


get_contrasts <- function(sampleData, pattern1, pattern2){
  first_contrasts <- unique(sampleData$group[grepl(pattern1, sampleData$group)])
  second_contrasts <- unique(sampleData$group[grepl(pattern2, sampleData$group)])
  first_contrasts <- setdiff(first_contrasts, second_contrasts)
  expand.grid(first_contrasts, second_contrasts)
}



get_results <- function(dds, c1, c2, outDir, prefix, andf){
  res <- results(dds, contrast=c('group', c1, c2), alpha=alpha, lfcThreshold=lfct, altHypothesis="greaterAbs")
  res$contrast  <-  paste0(c1, '_vs_', c2)
  res <- as.data.frame(res) %>% rownames_to_column("FLYBASE") %>% left_join(andf, by='FLYBASE')
  write.csv(res, file.path(outDir, paste0(prefix, c1, '_vs_', c2, '_l', lfct,'a', alpha, '_results.csv')),
            row.names=FALSE)
  return(res)
}


get_all_results <- function(sampleData, pattern1, pattern2, dds, outDir, prefix, annotations){
  contrasts <- get_contrasts(sampleData, pattern1, pattern2)
  f1 <- contrasts$Var1
  f2 <- contrasts$Var2
  andf <- read.csv(annotations, row.names=1)
  for(var in 1:length(f1)){
    print(paste0(as.character(f1[var]), ' vs. ', as.character(f2[var])))
    res <- get_results(dds,  as.character(f1[var]), as.character(f2[var]), outDir, prefix, andf)
  }
}
