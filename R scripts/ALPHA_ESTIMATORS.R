# ==============================================================================
#              ALPHA ESTIMATION METHODS FROM PIGA ET AL.
# ==============================================================================
## Derivative for NSB hyperprior
DerivativeNSB <- function(ns, alpha) {
  k <- length(ns)
  N <- sum(ns)
  s1 <- 0
  for (n in ns) {
    for (m in 1:n) {
      d <- m + alpha
      s1 <- s1 + 1/d
    }
  }
  s2 <- 0
  for (n in 1:N) {
    d <- n + k * alpha
    s2 <- s2 + k/d
  }
  q <- s1 - s2
  d2sbar <- k*k*psigamma(k*alpha+1, deriv=2)-psigamma(alpha+1, deriv=2)
  dsbar <- k*trigamma(k*alpha+1)-trigamma(alpha+1)
  return(q + d2sbar/dsbar)
}

## alpha with NSB hyperprior
alpha_starNSB <- function(ns, niter = 10) {
  alpha <- c(0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000)
  for (iter in 1:niter) {
    if (iter == 1) alphas <- alpha
    b <- alphas[1]
    q <- DerivativeNSB(ns, b)
    s <- -1
    if (q > 0) s <- 1
    for (b in alphas) {
      q <- DerivativeNSB(ns, b)
      if (q * s < 0) {
        if (iter == 1) {
          alphas <- sapply(1:10, function(n) b/10 + n*b/10)
          break
        } else {
          db <- alphas[2] - alphas[1]
          alphas <- sapply(1:10, function(n) b - db + db/10*n)
          break
        }
      }
    }
  }
  if (b == alpha[length(alpha)]) {
    if (DerivativeNSB(ns, b) < 0) {
      b <- alpha[1]
    }
  }
  return(list(b, DerivativeNSB(ns, b)))
}

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Function to calculate the derivative
Derivative <- function(vcounts, alpha) {
  k <- length(vcounts)
  N <- sum(vcounts)
  s1 <- 0
  for (n in vcounts) {
    for (m in 0:(n-1)) {
      d <- m + alpha
      s1 <- s1 + 1 / d
    }
  }
  s2 <- 0
  for (n in 0:(N-1)) {
    d <- n + k * alpha
    s2 <- s2 + k / d
  }
  q <- s1 - s2
  return(q)
}

# Function to estimate alpha_star
alpha_star <- function(vcounts, niter=10) {
  alpha <- c(0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000, 1000000)
  for (iter in 1:niter) {
    if (iter == 1) {
      alphas <- alpha
    }
    b <- alphas[1]
    q <- Derivative(vcounts, b)
    s <- -1
    if (q > 0) {
      s <- 1
    }
    for (b in alphas) {
      q <- Derivative(vcounts, b)
      if (q * s < 0) {
        if (iter == 1) {
          alphas <- sapply(0:10, function(n) b / 10 + n * b / 10)
          break
        } else {
          db <- alphas[2] - alphas[1]
          alphas <- sapply(0:10, function(n) b - db + db / 10 * n)
          break
        }
      }
    }
  }
  
  if (b == alpha[length(alpha)]) {
    if (Derivative(vcounts, b) < 0) {
      b <- alpha[1]
    }
  }
  return(list(b, Derivative(vcounts, b)))
}

# Function to calculate alpha_star for each column in a dataframe
calculate_alpha_star <- function(data, niter=10) {
  # Convert dataframe columns to a list of vectors
  vcounts_list <- lapply(data, function(col) as.numeric(col))
  
  # Calculate alpha_star for each vector of counts
  alpha_star_values <- sapply(vcounts_list, function(vcounts) {
    result <- alpha_star(vcounts, niter)
    return(result[[1]])
  })
  
  return(alpha_star_values)
}

# Example usage
# alpha_star_values <- calculate_alpha_star(kma_com)
# print(alpha_star_values)
# ==============================================================================
#                 ALPHA ESTIMATION USING MoM Method
# ==============================================================================
# Alpha calculation using the conventional MoM based on information from Wikipedia

# Define the function
compute_MoM_alpha <- function(data) {
  
  # Extract the name of the dataset
  count_set_name <- deparse(substitute(data))
  
  # Convert to matrix and remove the first column (assuming it is an index or non-numeric column)
  data_matrix <- as.matrix(data)
  
  # Calculate n (total counts per row)
  n <- rowSums(data, na.rm = TRUE)
  
  # Calculate mean proportions for each column
  mean_proportions <- colMeans(data / n, na.rm = TRUE)
  
  # Calculate variance proportions for each column
  variance_proportions <- apply(data / n, 2, function(x) var(x, na.rm = TRUE))
  
  # Handle division by zero
  variance_proportions[variance_proportions == 0] <- NA
  
  # Calculate alpha_j (for each component)
  alpha_j <- mean_proportions * (((mean_proportions * (1 - mean_proportions)) / variance_proportions) - 1)
  
  # Calculate alpha_0 (sum of alpha_j)
  alpha_0 <- sum(alpha_j, na.rm = TRUE)
  
  # Store results in the data frame
  results_df <- data.frame(
    Dataset = count_set_name,
    Min.Alpha = min(alpha_j, na.rm = TRUE),
    Max.Alpha = max(alpha_j, na.rm = TRUE),
    Mean.Alpha = mean(alpha_j, na.rm = TRUE),
    Median.Alpha = median(alpha_j, na.rm = TRUE),
    Alpha.Sum = alpha_0
  )
  
  # Remove names from alpha_j
  alphas <- unname(alpha_j)
  
  return(list(Summary = results_df, Alphas = alphas))
}

# ==============================================================================