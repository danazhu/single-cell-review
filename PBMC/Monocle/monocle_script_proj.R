#!/usr/bin/env Rscript

library(monocle)

tt = c()

DATASETS <- c("100", "1k")

for (ds in DATASETS) {
	dataset <- read.table(paste0('A_', ds, '.txt'))

	dataset <- as.matrix(dataset)

	start_time <- Sys.time()

	data <- newCellDataSet(dataset)
	data <- estimateSizeFactors(data)
	data <- reduceDimension(data, max_components = 3, num_dim = 10, reduction_method = 'tSNE', verbose = T, check_duplicates = F)
	data <- clusterCells(data, num_clusters = 10)

	end_time <- Sys.time()
	t = end_time - start_time
	print(t)

	write.csv(data$Cluster, sprintf(paste0("monocle_output_", ds, ".txt")));
}
