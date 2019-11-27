#!/usr/bin/env Rscript
setwd("~/Desktop/GroupProject/single-cell-review/PBMC/pcaReduce")
library("pcaReduce")

# Source ARI function
source("../adjusted_rand_index.R")
# Read in ground truth files
library(rhdf5)
y_100 = h5read("../A_samples.mat", "y_100") 
y_1k = h5read("../A_samples.mat", "y_1k")

tt = c()

DATASETS <- c("100", "1k")
TRIALS <- c(1:5)

# in case ari and time output file exists already
file.remove("pca_time_output.txt")
file.remove("pca_ari_output.txt")

for (ds in DATASETS) {
  dataset <- read.table(paste0('A_', ds, '.txt'))
  dataset <- as.matrix(dataset)
  
  write.table(paste0("Time values for ", ds, " sample size:"), "pca_time_output.txt", append = TRUE, sep = "\n", col.names = FALSE, row.names = FALSE)  
  
  for (num in TRIALS) {
    D <- log2(dataset + 1) 
    Input <- t(D) 
    start_time <- Sys.time()
    Output_S <- PCAreduce(Input, nbt=1, q=9, method='S') 
    end_time <- Sys.time()
    
    t = end_time - start_time 
    
    # calculate ari values
    c <- matrix(Output_S[[1]], nrow = length(Output_S[[1]]), ncol = 1)
    if (ds == "100") {
      current_ARI <- compute_ARI(c, y_100)
    }
    else if (ds == "1k") {
      current_ARI <- compute_ARI(c, y_1k)
    }
    
    # record time spent
    write.table(end_time - start_time, "pca_time_output.txt", append = TRUE, sep="\n", col.names = FALSE, row.names = FALSE)
    write.table(current_ARI, "pca_ari_output.txt", append = TRUE, sep="\n", col.names = FALSE, row.names = FALSE)
    
    write.csv(Output_S[1], sprintf(paste0("pca_output_", ds, "_", num, ".txt")));
  }
}
