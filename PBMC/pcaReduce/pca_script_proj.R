#!/usr/bin/env Rscript
library("pcaReduce")

tt = c()

DATASETS <- c("100", "1k")

for (ds in DATASETS) {
  dataset <- read.table(paste0('A_', ds, '.txt'))
  dataset <- as.matrix(dataset)
  
  start_time <- Sys.time()
  
  D <- log2(dataset + 1) 
  Input <- t(D) 
  start_time <- Sys.time()
  Output_S <- PCAreduce(Input, nbt=1, q=9, method='S')
  end_time <- Sys.time()
  
  t = end_time - start_time 
  
  print(t)
  
  write.csv(Output_S[1], sprintf(paste0("pca_output_", ds, ".txt")));
}

