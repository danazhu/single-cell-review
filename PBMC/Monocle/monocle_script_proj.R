#!/usr/bin/env Rscript

library(monocle)

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
file.remove("monocle_time_output.txt")
file.remove("monocle_ari_output.txt")

for (ds in DATASETS) {
  dataset <- read.table(paste0('A_', ds, '.txt'))
  dataset <- as.matrix(dataset)
  
  write.table(paste0("time values for ", ds, " sample size"), "monocle_time_output.txt", append = TRUE, sep = "\n", col.names = FALSE, row.names = FALSE)
  
  for (num in TRIALS) {
    start_time <- Sys.time()
    
    data <- newCellDataSet(dataset)
    data <- estimateSizeFactors(data)
    data <- reduceDimension(data, max_components = 3, num_dim = 10, reduction_method = 'tSNE', verbose = T, check_duplicates = F)
    data <- clusterCells(data, num_clusters = 10)
    
    end_time <- Sys.time()
    t = end_time - start_time
    #print(t)
    print(end_time - start_time)
    
    # calculate ari values
    c <- matrix(as.numeric(data$Cluster), nrow = length(data$Cluster), ncol = 1)
    if (ds == "100") {
      current_ARI <- compute_ARI(c, y_100)
    }
    else if (ds == "1k") {
      current_ARI <- compute_ARI(c, y_1k)
    }
    
    # record time spent
    write.table(end_time - start_time, "monocle_time_output.txt", append = TRUE, sep="\n", col.names = FALSE, row.names = FALSE)
    write.table(current_ARI, "monocle_ari_output.txt", append = TRUE, sep="\n", col.names = FALSE, row.names = FALSE)
    
    write.csv(data$Cluster, sprintf(paste0("monocle_output_", ds, "_", num, ".txt")));
  }
}