# Argument errors (optional)
p1 <- c(0, 1, 3, 4, 5, 2, 5, 5, 5)
p2 <- c(0, 1, 3, 4, 5, 2, 5, 3, 4)
RI <- compute_ARI(p1, p2)

# Manually compute adjusted RI
compute_ARI <- function(p1, p2) {
  # Find unique partition lengths (number of clusters)
  N <- length(p1) # length of partition set
  N1 <- max(unique(p1))
  N2 <- max(unique(p2))
  
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
