# FUNCTIONS FOR ALL ENTROPY ESTIMATORS USED IN THE THESIS PROJECT
# ==============================================================================
source("ALPHA_ESTIMATORS.R")
# ==============================================================================
# CHAO-SHEN ESTIMATOR -- BASED ON the entropy::entropy.ChaoShen() FUNCTION
# ==============================================================================
chao_shen_entropy <- function(y, unit = c("log", "log2", "log10")) {
  unit = match.arg(unit)
  yx = y[y > 0] # Filter out zero counts
  N = sum(yx) # Total count
  p_MLE = yx / N # Proportions
  f1 = sum(yx == 1) # Singletons
  if (f1 == N) {
    f1 = N - 1
  }
  C = 1 - f1 / N # Coverage
  C_hat = C * p_MLE # Adjusted proportions
  la = (1 - (1 - C_hat)^N) # Probability adjustment
  H_CS = -sum(C_hat * log(C_hat) / la) # Entropy calculation
  if (unit == "log2") {
    H_CS = H_CS / log(2)
  }
  if (unit == "log10") {
    H_CS = H_CS / log(10)
  }
  return(H_CS)
}
# ==============================================================================
# HAUSSER-STRIMMER ESTIMATOR
# ==============================================================================
hausser_strimmer_entropy <- function(counts) {
  # Total count of all elements
  N <- sum(counts)
  
  # Total number of unique elements: taxa or species categories
  K <- ncol(counts)
  
  # Maximum likelihood estimates for each category
  theta_ML <- t(colSums(counts)) / N
  
  # Shrinkage target
  t_k <- 1 / K
  
  # Calculate lambda
  lambda <- (1 - sum((colSums(counts) / N)^2)) / ((N - 1) * sum((t_k - colSums(counts) / N)^2))
  lambda <- pmax(0, pmin(1, lambda))  # Ensuring lambda is within [0,1]
  
  # Shrinkage estimator for each category
  theta_Shrink <- lambda * t_k + (1 - lambda) * theta_ML
  
  # Compute Hausser-Strimmer estimator
  H_HS <- -sum(theta_Shrink * log(theta_Shrink))
  
  return(H_HS)
}

# OR
# Shrinkage intensity
lambda_hs <- function(multi, t = 1/100) {
  K <- length(multi)
  N <- sum(multi)
  return((1 - sum((multi/N)^2)) / ((N - 1) * sum((t - multi/N)^2)))
}

hs_entropy <- function(counts) {
  multi <- colSums(counts)
  K <- length(multi)
  N <- sum(multi)
  t <- 1 / K # shrinkage target
  l <- lambda_hs(multi, t)
  l <- max(0, min(1, l))
  rhos <- l*t + (1 - l) * multi/N # multi/N is the maximum likelihood estimate
  return(list(entropy(rhos)))
}
# ==============================================================================
# WOLPERT-WOLF ESTIMATOR
# ==============================================================================
# Entropy estimated from counts using wolpert-wolf estimator and MoM alpha
ww_entropy <- function(counts) {
  ns <- colSums(counts)
  alpha <- compute_MoM_alpha(counts)[[2]] # compute_MoM_alpha function from ALPHA_ESTIMATORS.R
  alpha_tilde <- ns + alpha
  A_tilde <- sum(alpha_tilde)
  ww_ent <- digamma(A_tilde + 1) - sum((alpha_tilde / A_tilde) * digamma(alpha_tilde + 1))
  return(ww_ent)
}

# For NSB and Piga methods
aww_entropy <- function(ns, alpha) {
  sapply(alpha, function(a) {
    alpha_tilde <- ns + a
    A_tilde <- sum(alpha_tilde)
    ww_ent <- digamma(A_tilde + 1) - sum((alpha_tilde / A_tilde) * digamma(alpha_tilde + 1))
    return(ww_ent)
  })
}
# ==============================================================================
# NSB ESTIMATOR
# ==============================================================================
# NSB hyperprior using the trigamma function
hyperprior_nsb <- function(alpha, K) {
  numerator <- K * trigamma(K * alpha + 1) - trigamma(alpha + 1)
  denominator <- log(K)
  p_alpha <- numerator / denominator
  return(p_alpha)
}

# Aggregated entropy for all category counts
nsb_entropy <- function(counts) {
  ns <- colSums(counts)
  K <- length(ns)
  # Numerical integration over alpha value from 0 to infinity
  result <- integrate(
    function(alpha) {
      ent <- aww_entropy(ns, alpha)
      hyp <- hyperprior_nsb(alpha, K)
      ent * hyp
    },
    lower = 0, upper = Inf, subdivisions = 1000, rel.tol = .Machine$double.eps^0.25
  )
  if (result$message != "OK") {
    warning("Integration may not be accurate:", result$message)
  }
  return(result$value)
}

# NSB single to be used only with the apply function when entropy for each sample is required
# Also can be used with colSums vector of data
nsb_entropy_single <- function(counts) {
  ns <- counts
  K <- length(ns)
  # Numerical integration over alpha value from 0 to infinity
  result <- integrate(
    function(alpha) {
      ent <- aww_entropy(ns, alpha)
      hyp <- hyperprior_nsb(alpha, K)
      ent * hyp
    },
    lower = 0, upper = Inf, subdivisions = 1000, rel.tol = .Machine$double.eps^0.25
  )
  if (result$message != "OK") {
    warning("Integration may not be accurate:", result$message)
  }
  return(result$value)
}

# nsb_entropy_single can be used iteratively for each sample by using:
# apply(data, 1, function(row) nsb_entropy_single(as.numeric(row)))
# nsb_entropy(data)
# ==============================================================================
# PIGA ESTIMATOR
# ==============================================================================
piga_entropy <- function(counts) {

  # Estimate values for ns, K and N from the dataset
  ns <- colSums(counts)
  K <- length(ns)
  N <- sum(counts)

  # Estimate alpha using alpha Star function
  estimated_alpha <- alpha_star(ns, niter = 10)
  alpha <- as.numeric(estimated_alpha[[1]])

  # Calculate entropy using moment given in Piga's paper
  first_term <- digamma(N + K * alpha + 1)
  second_term <- sum(((ns + alpha) / (N + K * alpha)) * digamma(ns + alpha + 1))
  piga_entropy_estimate <-  first_term - second_term

  return(piga_entropy_estimate)
}

# Piga entropy for a single sample only to be used with the apply function
# Also can be used with colSums vector of data
piga_entropy_single <- function(counts) {
  
  # # Estimate values for ns, K and N from the dataset -- Do not do a colSums
  ns <- counts
  K <- length(ns)
  N <- sum(counts)
  
  # Estimate alpha using alpha Star function
  estimated_alpha <- alpha_star(ns, niter = 10)
  alpha <- as.numeric(estimated_alpha[[1]])
  
  # Calculate entropy using moment given in Piga's paper
  first_term <- digamma(N + K * alpha + 1)
  second_term <- sum(((ns + alpha) / (N + K * alpha)) * digamma(ns + alpha + 1))
  piga_entropy_estimate <-  first_term - second_term
  
  return(piga_entropy_estimate)
}

# How to use to get entropy for each sample
# apply(data, 1, function(row) piga_entropy_single(as.numeric(row)))
# ==============================================================================
# Function for estimating entropy for aggregated counts for an abundance dataset
estimate_aggregated_entropies <- function(data) {
  
  # Shannon's entropy
  H_shannon <- diversity(colSums(data), index = "shannon")
  
  # Chao-Shen entropy
  H_cs <- entropy.ChaoShen(colSums(data))
  
  
  # Hausser-Strimmer entropy
  H_hs <- entropy.shrink(colSums(data), verbose = FALSE)
  
  
  # Wolpert-Wolf entropy
  H_ww <- ww_entropy(data)
  
  # NSB entropy
  H_nsb <- nsb_entropy_single(colSums(data))
  
  # Piga -- piga_entropy_single gives a higher value with colSums than piga_entropy
  H_piga <- piga_entropy_single(colSums(data))
  
  # Create a dataframe of all entropy estimate
  entropies_df <- data.frame(
    "Naive.entropy" = H_shannon,
    "Chao.Shen" = H_cs,
    "Hausser.Strimmer" = H_hs,
    "Wolpert.Wolf" = H_ww,
    "NSB" = H_nsb,
    "Piga" = H_piga)
  
  return(entropies_df)
}
# ==============================================================================
# Function for estimating entropy for each sample
estimate_individual_entropies <- function(data) {
  
  # Shannon's entropy
  H_shannon <- H_shannon <- diversity(data, index="shannon")
  
  # Chao-Shen entropy
  H_cs <- apply(data, 1, function(row) entropy.ChaoShen(as.numeric(row)))
  
  # Hausser-Strimmer entropy
  H_hs <- apply(data, 1, function(row) entropy.shrink(as.numeric(row), verbose=FALSE))
  
  # Wolpert-Wolf entropy
  alpha_values <- compute_MoM_alpha(data)[[2]]
  ww_alpha <- median(alpha_values, na.rm = TRUE) # Median alpha used to reduce impact of extreme values
  #ww_alpha <- mean(alpha_values, na.rm = TRUE) 
  if (is.na(ww_alpha)) ww_alpha <- 0 
  H_ww <- apply(data, 1, function(row) aww_entropy(as.numeric(row), ww_alpha))
  
  # NSG entropy
  H_nsb <- apply(data, 1, function(row) nsb_entropy_single(as.numeric(row)))
  
  # Piga entropy
  H_piga <- apply(data, 1, function(row) piga_entropy_single(as.numeric(row)))
  
  # Create entropy dataframe
  entropies_df <- data.frame(
    "Naive.entropy" = H_shannon,
    "Chao.Shen" = H_cs,
    "Hausser.Strimmer" = H_hs,
    "Wolpert.Wolf" = H_ww,
    "NSB" = H_nsb,
    "Piga" = H_piga)
  
  return(entropies_df)
}
# ==============================================================================
