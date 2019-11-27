# Read in ground truth tables of clusters
library(rhdf5)
y_100 = h5read("A_samples.mat", "y_100") 
y_1k = h5read("A_samples.mat", "y_1k")

# Sample: Test Monocle 100 result
p1 <- read.delim("single-cell-review/PBMC/Monocle/monocle_output_100_1.txt", header = TRUE, sep = ",")
p1 <- p1$x
p1 <- matrix(p1, nrow = length(p1), ncol = 1)

RI_monocle_100 <- compute_ARI(p1, y_100)

# Sample: Test pcaReduce 100 result
pca1 <- read.delim("single-cell-review/PBMC/pcaReduce/pca_output_100_1.txt", header = TRUE, sep = ",")
pca1 <- pca1$Cl_mat
pca1 <- matrix(pca1, nrow = length(pca1), ncol = 1)

RI_pca_100 <- compute_ARI(pca1, y_100)

# Sample: Test Monocle 1k result
p2 <- read.delim("single-cell-review/PBMC/Monocle/monocle_output_1k_4.txt", header = TRUE, sep = ",")
p2 <- p2$x
p2 <- matrix(p2, nrow = length(p2), ncol = 1)
RI_monocle_1k <- compute_ARI(p2, y_1k)

# Sample: Test pcaReduce 1k result
pca2 <- read.delim("single-cell-review/PBMC/pcaReduce/pca_output_1k_1.txt", header = TRUE, sep = ",")
pca2 <- pca2$Cl_mat
pca2 <- matrix(pca2, nrow = length(pca2), ncol = 1)

RI_pca_1k <- compute_ARI(pca2, y_1k)

# Manually compute adjusted RI
compute_ARI <- function(p1, p2) {
  # Find unique partition lengths (number of clusters)
  N <- length(p1) # length of partition set
  # Mimic the unique() function in MATLAB
  # Find set of unique, sorted values
  p1_unique_sort <- sort(unique(p1))
  p2_unique_sort <- sort(unique(p2))
  # Convert p1, p2 values into indexes of unique sort
  p1_indices <- apply(p1, c(1, 2), function(elt) result <- match(elt, p1_unique_sort))
  p2_indices <- apply(p2, c(1, 2), function(elt) result <- match(elt, p2_unique_sort))
  
  p1 <- p1_indices
  p2 <- p2_indices
  N1 <- max(p1)
  N2 <- max(p2)
  
  n <- matrix(-1, nrow = N1, ncol = N2)
  # Create n)i,j the match matrix
  for (i in 1:N1) {
    for (j in 1:N2) {
      G1 <- which(p1 == i)
      G2 <- which(p2 == j)
      
      # assign how many p1_i and p2_j have in common
      n[i, j] <- length(intersect(G1, G2))
      
    }
  }
  
  # Calculate adjusted rand index
  ssm <- 0
  sm1 <- 0
  sm2 <- 0
  
  # Calculate ssm
  for (i in 1:N1) {
    for (j in 1:N2) {
      ssm <- ssm + nchoose2(n[i, j])
    }
  }
  
  # Calculate sm1 (across columns, per row horizontally)
  a_i <- apply(n, 1, sum)
  for (i in 1:N1) {
    sm1 <- sm1 + nchoose2(a_i[i])
  }
  
  # Calculate sm2 (across rows, per column vertically)
  b_j <- apply(n, 2, sum)
  for (i in 1:N2) {
    sm2 <- sm2 + nchoose2(b_j[i])
  }
  
  NN <- ssm - sm1 * sm2 / nchoose2(N)
  DD <- (sm1 + sm2) / 2 - sm1 * sm2 / nchoose2(N)
  RI <- NN/DD
}

# Calculate positive combination (n 2)
nchoose2 <- function(n) {
  if (n > 1) {
    result <- choose(n, 2)
  } else {
    result <- 0
  }
}
