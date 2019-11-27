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
