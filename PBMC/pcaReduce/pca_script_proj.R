#!/usr/bin/env Rscript
library("pcaReduce")

tt = c()

DATASETS <- c("100", "1k")
TRIALS <- c(1:5)

for (ds in DATASETS) {
  dataset <- read.table(paste0('A_', ds, '.txt'))
  dataset <- as.matrix(dataset)
  
  write.table(paste0("ari values for ", ds, " sample size"), "pca_ari_time_output.txt", append = TRUE, sep = "\n", col.names = FALSE, row.names = FALSE)
  
  for (num in TRIALS) {
    D <- log2(dataset + 1) 
    Input <- t(D) 
    start_time <- Sys.time()
    Output_S <- PCAreduce(Input, nbt=1, q=9, method='S') 
    end_time <- Sys.time()
    
    t = end_time - start_time 
    
    # TODO: calculate ari values, mean run time between the 5 samples
    
    # record time spent
    write.table(end_time - start_time, "pca_ari_time_output.txt", append = TRUE, sep="\n", col.names = FALSE, row.names = FALSE)
    
    write.csv(Output_S[1], sprintf(paste0("pca_output_", ds, "_", num, ".txt")));
  }
}
