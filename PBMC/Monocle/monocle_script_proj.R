#!/usr/bin/env Rscript

library(monocle)

tt = c()

DATASETS <- c("100", "1k")
TRIALS <- c(1:5)

# in case ari and time output file exists already
file.remove("monocle_ari_time_output.txt")

for (ds in DATASETS) {
	dataset <- read.table(paste0('A_', ds, '.txt'))
	dataset <- as.matrix(dataset)

	write.table(paste0("ari values for ", ds, " sample size"), "monocle_ari_time_output.txt", append = TRUE, sep = "\n", col.names = FALSE, row.names = FALSE)

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
		
		# TODO: calculate ari values

		# record time spent
		write.table(end_time - start_time, "monocle_ari_time_output.txt", append = TRUE, sep="\n", col.names = FALSE, row.names = FALSE)

		write.csv(data$Cluster, sprintf(paste0("monocle_output_", ds, "_", num, ".txt")));
	}
}
